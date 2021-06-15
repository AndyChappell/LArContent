/**
 *  @file   larpandoracontent/LArMonitoring/ShowerCaloMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the shower calo algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/ShowerCaloMonitoringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerCaloMonitoringAlgorithm::ShowerCaloMonitoringAlgorithm() : m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerCaloMonitoringAlgorithm::~ShowerCaloMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerCaloMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    this->MakeSelection(pMCParticleList, pCaloHitList, targetMCParticleToHitsMap);

    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (auto [pMCParticle, caloHits] : targetMCParticleToHitsMap)
    {
        if (!(pMCParticle->GetParticleId() == PHOTON || std::abs(pMCParticle->GetParticleId()) == E_MINUS))
        {
            continue;
        }
        else
        {
            const CartesianVector &vertex(pMCParticle->GetVertex());
            const CartesianVector &endpoint(pMCParticle->GetEndpoint());
            const CartesianVector startU(vertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(vertex.GetY(), vertex.GetZ())));
            const CartesianVector startV(vertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(vertex.GetY(), vertex.GetZ())));
            const CartesianVector startW(vertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(vertex.GetY(), vertex.GetZ())));
            const CartesianVector finishU(endpoint.GetX(), 0.f, static_cast<float>(transform->YZtoU(endpoint.GetY(), endpoint.GetZ())));
            const CartesianVector finishV(endpoint.GetX(), 0.f, static_cast<float>(transform->YZtoV(endpoint.GetY(), endpoint.GetZ())));
            const CartesianVector finishW(endpoint.GetX(), 0.f, static_cast<float>(transform->YZtoW(endpoint.GetY(), endpoint.GetZ())));
            const CartesianVector &axisU{(finishU - startU).GetUnitVector()};
            const CartesianVector &axisV{(finishV - startV).GetUnitVector()};
            const CartesianVector &axisW{(finishW - startW).GetUnitVector()};

            unsigned int nHitsU{0}, nHitsV{0}, nHitsW{0};
            CaloHitList wHits;
            for (const CaloHit *pCaloHit : caloHits)
            {
                const CartesianVector &position{pCaloHit->GetPositionVector()};
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                    {
                        ++nHitsU;
                        const CartesianVector &rel{position - startU};
                        const float dot{rel.GetDotProduct(axisU)};
                        CartesianVector proj(axisU);
                        proj *= dot;
                        proj += startU;
                        break;
                    }
                    case TPC_VIEW_V:
                    {
                        ++nHitsV;
                        const CartesianVector &rel{position - startV};
                        const float dot{rel.GetDotProduct(axisV)};
                        CartesianVector proj(axisV);
                        proj *= dot;
                        proj += startV;
                        break;
                    }
                    case TPC_VIEW_W:
                    {
                        wHits.emplace_back(pCaloHit);
                        ++nHitsW;
                        const CartesianVector &rel{position - startW};
                        const float dot{rel.GetDotProduct(axisW)};
                        CartesianVector proj(axisW);
                        proj *= dot;
                        proj += startW;
                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &proj, "Hit", BLACK, 1);
                        break;
                    }
                    default:
                        break;
                }
            }
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "W", RED));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerCaloMonitoringAlgorithm::MakeSelection(
    const MCParticleList *pMCList, const CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcMap) const
{
    // Default reconstructability criteria are very liberal to allow for unfolded hierarchy
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 15;
    parameters.m_minHitsForGoodView = 5;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_minHitSharingFraction = 0;
    parameters.m_foldBackHierarchy = false;

    LArMCParticleHelper::SelectReconstructableMCParticles(pMCList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerCaloMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
