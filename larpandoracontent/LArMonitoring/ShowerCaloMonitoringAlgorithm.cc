/**
 *  @file   larpandoracontent/LArMonitoring/ShowerCaloMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the shower calo algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/ShowerCaloMonitoringAlgorithm.h"

#include <memory>

using namespace pandora;

namespace lar_content
{

ShowerCaloMonitoringAlgorithm::ShowerCaloMonitoringAlgorithm() : m_visualize(false), m_writeTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerCaloMonitoringAlgorithm::~ShowerCaloMonitoringAlgorithm()
{
    if (m_writeTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "ShowerCaloMonitoringAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerCaloMonitoringAlgorithm::Run()
{
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    this->MakeSelection(pMCParticleList, pCaloHitList, targetMCParticleToHitsMap);

    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (auto [pMCParticle, caloHits] : targetMCParticleToHitsMap)
    {
        const int pdg{std::abs(pMCParticle->GetParticleId())};
        const int tier{LArMCParticleHelper::GetHierarchyTier(pMCParticle)};
        if (!(pdg == PHOTON || pdg == E_MINUS))
        {
            continue;
        }
        else
        {
            std::cout << "E(MC): " << pMCParticle->GetEnergy() << std::endl;
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
            ProjCaloHitPtrList projCaloHitListU, projCaloHitListV, projCaloHitListW;
            for (const CaloHit *pCaloHit : caloHits)
            {
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                    {
                        ++nHitsU;
                        this->ProjectHit(pCaloHit, startU, axisU, projCaloHitListU);
                        break;
                    }
                    case TPC_VIEW_V:
                    {
                        ++nHitsV;
                        this->ProjectHit(pCaloHit, startV, axisV, projCaloHitListV);
                        break;
                    }
                    case TPC_VIEW_W:
                    {
                        ++nHitsW;
                        this->ProjectHit(pCaloHit, startW, axisW, projCaloHitListW);
                        break;
                    }
                    default:
                        break;
                }
            }
            std::sort(projCaloHitListW.begin(), projCaloHitListW.end(),
                [](const ProjCaloHitPtr pLhs, const ProjCaloHitPtr pRhs) { return pLhs->first < pRhs->first; });
            if (m_visualize && !projCaloHitListW.empty())
            {
                const CartesianVector &start(projCaloHitListW.front()->second->GetPositionVector());
                const CartesianVector &finish(projCaloHitListW.back()->second->GetPositionVector());
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &finish, "Projection", BLUE, 1, 1));

                CaloHitList wHits;
                for (const ProjCaloHitPtr pProjCaloHit : projCaloHitListW)
                {
                    std::cout << "Pos: " << pProjCaloHit->first << " ADC: " << pProjCaloHit->second->GetElectromagneticEnergy() << std::endl;
                    wHits.emplace_back(pProjCaloHit->second);
                }
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "W", RED));
            }
            if (m_writeTree)
            {
                float mip{0.f}, em{0.f}, had{0.f};
                FloatVector projPositionVectorU, adcVectorU;
                for (const ProjCaloHitPtr pProjCaloHit : projCaloHitListU)
                {
                    projPositionVectorU.emplace_back(pProjCaloHit->first);
                    adcVectorU.emplace_back(pProjCaloHit->second->GetElectromagneticEnergy());
                    mip += pProjCaloHit->second->GetMipEquivalentEnergy();
                    em += pProjCaloHit->second->GetElectromagneticEnergy();
                    had += pProjCaloHit->second->GetHadronicEnergy();
                }
                FloatVector projPositionVectorV, adcVectorV;
                for (const ProjCaloHitPtr pProjCaloHit : projCaloHitListV)
                {
                    projPositionVectorV.emplace_back(pProjCaloHit->first);
                    adcVectorV.emplace_back(pProjCaloHit->second->GetElectromagneticEnergy());
                    mip += pProjCaloHit->second->GetMipEquivalentEnergy();
                    em += pProjCaloHit->second->GetElectromagneticEnergy();
                    had += pProjCaloHit->second->GetHadronicEnergy();
                }
                FloatVector projPositionVectorW, adcVectorW;
                for (const ProjCaloHitPtr pProjCaloHit : projCaloHitListW)
                {
                    projPositionVectorW.emplace_back(pProjCaloHit->first);
                    adcVectorW.emplace_back(pProjCaloHit->second->GetElectromagneticEnergy());
                    mip += pProjCaloHit->second->GetMipEquivalentEnergy();
                    em += pProjCaloHit->second->GetElectromagneticEnergy();
                    had += pProjCaloHit->second->GetHadronicEnergy();
                }

                std::cout << "#EM: " << em << " #Had: " << had << " #MIP: " << mip << std::endl;

                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pdg", pdg));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "tier", tier));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "projectedPositionsU", &projPositionVectorU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energiesU", &adcVectorU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "projectedPositionsV", &projPositionVectorV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energiesV", &adcVectorV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "projectedPositionsW", &projPositionVectorW));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energiesW", &adcVectorW));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }
        }
        if (m_visualize)
        {
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
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

void ShowerCaloMonitoringAlgorithm::ProjectHit(
    const CaloHit *pCaloHit, const CartesianVector &origin, const CartesianVector &axis, ProjCaloHitPtrList &projCaloHitList)
{
    const CartesianVector &rel{pCaloHit->GetPositionVector() - origin};
    const float dot{rel.GetDotProduct(axis)};
    projCaloHitList.emplace_back(ProjCaloHitPtr(new ProjCaloHit(dot, pCaloHit)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerCaloMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    }


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
