/**
 *  @file   larpandoracontent/LArMonitoring/PfoValidationAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/PfoValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoValidationAlgorithm::PfoValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoValidationAlgorithm::~PfoValidationAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "pfo_tree", "pfo_tree.root", "RECREATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "PfoValidationAlgorithm: Unable to write pfo_tree to file" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoValidationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    for (std::string pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        std::map<const MCParticle *, float> hitMap;
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            try
            {
                const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (!pMCParticle)
                    continue;
                const float adc{pCaloHit->GetInputEnergy()};
                hitMap[pMCParticle] += adc;
            }
            catch (const StatusCodeException &)
            {
            }
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            CaloHitList caloHitListU, caloHitListV, caloHitListW, caloHitList;
            LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_U, caloHitListU);
            LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_V, caloHitListV);
            LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_W, caloHitListW);
            LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_U, caloHitList);
            LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_V, caloHitList);
            LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_W, caloHitList);
            int nGoodViews{0};
            nGoodViews += caloHitListU.size() >= 5 ? 1 : 0;
            nGoodViews += caloHitListV.size() >= 5 ? 1 : 0;
            nGoodViews += caloHitListW.size() >= 5 ? 1 : 0;
            const int nTotalHits{static_cast<int>(caloHitList.size())};
            if ((nGoodViews < 3) || (nTotalHits < 15))
                continue;

            float showerContribution{0.f}, trackContribution{0.f}, totalContribution{0.f}, neutronContribution{0.f}, recoAdc{0.f};
            std::map<const MCParticle *, float> contributionMap;
            for (const CaloHit *pCaloHit : caloHitList)
            {
                try
                {
                    const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                    if (!pMCParticle)
                        continue;
                    const int truePdg{static_cast<int>(std::abs(pMCParticle->GetParticleId()))};
                    const float adc{pCaloHit->GetInputEnergy()};
                    const MCParticle *pParent{pMCParticle};
                    bool isDownstreamNeutron{false};
                    contributionMap[pMCParticle] += adc;
                    while (!pParent->GetParentList().empty())
                    {
                        pParent = pParent->GetParentList().front();
                        if (std::abs(pParent->GetParticleId()) == NEUTRON)
                        {
                            neutronContribution += adc;
                            isDownstreamNeutron = true;
                            break;
                        }
                    }
                    recoAdc += adc;
                    totalContribution += adc;
                    if (isDownstreamNeutron)
                        continue;
                    if (truePdg == E_MINUS || truePdg == PHOTON)
                        showerContribution += adc;
                    else
                        trackContribution += adc;
                }
                catch (const StatusCodeException &)
                {
                }
            }
            if (totalContribution <= std::numeric_limits<float>::epsilon())
                continue;
            const bool isDownstreamNeutron{(neutronContribution > trackContribution) && (neutronContribution > showerContribution) ? true : false};
            if (isDownstreamNeutron)
                continue;
            const MCParticle *pBestMC{nullptr};
            float maxContribution{0.f};
            for (const auto &[ pMC, contribution ] : contributionMap)
            {
                if (contribution > maxContribution)
                {
                    maxContribution = contribution;
                    pBestMC = pMC;
                }
            }
            const int bestPdg{std::abs(pBestMC->GetParticleId())};
            const float trueEnergy{pBestMC->GetEnergy()};
            const int isRecoTrack{pPfo->GetParticleId() == MU_MINUS ? 1 : 0};
            const int isTrueTrack{bestPdg == E_MINUS || bestPdg == PHOTON ? 0 : 1};
            const float purity{maxContribution / totalContribution};
            const float completeness{maxContribution / hitMap[pBestMC]};
         
            std::cout << "PDG " << bestPdg << " Et " << trueEnergy << " ADCt " << contributionMap[pBestMC] << " ADCr " << recoAdc <<
                " P: " << purity << " C: " << completeness << " Track_t " << isTrueTrack << " Track_r " << isRecoTrack << std::endl;
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "true_pdg", bestPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "true_energy", trueEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "true_adc", contributionMap[pBestMC]));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "reco_adc", recoAdc));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "purity", purity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "completeness", completeness));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "is_true_track", isTrueTrack));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pfo_tree", "is_reco_track", isRecoTrack));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), "pfo_tree"));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
