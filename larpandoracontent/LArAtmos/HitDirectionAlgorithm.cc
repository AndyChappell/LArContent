/**
 *  @file   larpandoracontent/LArMonitoring/HitDirectionAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArAtmos/HitDirectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

HitDirectionAlgorithm::HitDirectionAlgorithm() :
    m_writeFile{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitDirectionAlgorithm::~HitDirectionAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitDirectionAlgorithm::Run()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
    if (!primaries.empty())
    {
        for (const MCParticle *primary : primaries)
        {
            const MCParticleList &parents{primary->GetParentList()};
            if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
            {
                pTrueNeutrino = parents.front();
                break;
            }
        }
    }

    if (pTrueNeutrino)
    {
        double uAdcTotal{0.}, vAdcTotal{0.}, wAdcTotal{0.};
        double uXAvg{0.}, uZAvg{0.}, vXAvg{0.}, vZAvg{0.}, wXAvg{0.}, wZAvg{0.};
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            const HitType view{pCaloHit->GetHitType()};
            const CartesianVector position{pCaloHit->GetPositionVector()};
            const double adc{pCaloHit->GetInputEnergy()};
            switch (view)
            {
                case TPC_VIEW_U:
                    uAdcTotal += adc;
                    uXAvg += adc * position.GetX();
                    uZAvg += adc * position.GetZ();
                    break;
                case TPC_VIEW_V:
                    vAdcTotal += adc;
                    vXAvg += adc * position.GetX();
                    vZAvg += adc * position.GetZ();
                    break;
                case TPC_VIEW_W:
                    wAdcTotal += adc;
                    wXAvg += adc * position.GetX();
                    wZAvg += adc * position.GetZ();
                    break;
                default:
                    break;
            }
        }
        if (uAdcTotal > std::numeric_limits<double>::epsilon())
        {
            uXAvg /= uAdcTotal;
            uZAvg /= uAdcTotal;
        }
        if (vAdcTotal > std::numeric_limits<double>::epsilon())
        {
            vXAvg /= vAdcTotal;
            vZAvg /= vAdcTotal;
        }
        if (wAdcTotal > std::numeric_limits<double>::epsilon())
        {
            wXAvg /= wAdcTotal;
            wZAvg /= wAdcTotal;
        }

        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
        const CartesianVector &trueMomentum{pTrueNeutrino->GetMomentum()};
        const CartesianVector &momOffset{trueMomentum.GetUnitVector() + trueVertex};
        const CartesianVector vu(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ())));
        const CartesianVector vv(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ())));
        const CartesianVector vw(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ())));
        const CartesianVector momu(momOffset.GetX(), 0.f, static_cast<float>(transform->YZtoU(momOffset.GetY(), momOffset.GetZ())));
        const CartesianVector momv(momOffset.GetX(), 0.f, static_cast<float>(transform->YZtoV(momOffset.GetY(), momOffset.GetZ())));
        const CartesianVector momw(momOffset.GetX(), 0.f, static_cast<float>(transform->YZtoW(momOffset.GetY(), momOffset.GetZ())));
        const CartesianVector pu{momu - vu};
        const CartesianVector pv{momv - vv};
        const CartesianVector pw{momw - vw};

        const double u_mom_x{uXAvg - trueVertex.GetX()}, u_mom_z{uZAvg - vu.GetZ()};
        const double v_mom_x{vXAvg - trueVertex.GetX()}, v_mom_z{vZAvg - vv.GetZ()};
        const double w_mom_x{wXAvg - trueVertex.GetX()}, w_mom_z{wZAvg - vw.GetZ()};

        std::cout << "True Direction 3D: " << trueMomentum.GetUnitVector() << std::endl;
        const double norm{std::sqrt(w_mom_x * w_mom_x + w_mom_z * w_mom_z)};
        std::cout << "True Direction W: " << (w_mom_x / norm) << ", " << (w_mom_z / norm) << std::endl;

        const int truePdg{pTrueNeutrino->GetParticleId()};
        if (m_writeFile)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_pdg", truePdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_mom_x", trueMomentum.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_mom_y", trueMomentum.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_mom_z", trueMomentum.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_mom_u", pu.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_mom_v", pv.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_mom_w", pw.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "u_mom_x", u_mom_x));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "u_mom_z", u_mom_z));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "v_mom_x", v_mom_x));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "v_mom_z", v_mom_z));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "w_mom_x", w_mom_x));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "w_mom_z", w_mom_z));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitDirectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
