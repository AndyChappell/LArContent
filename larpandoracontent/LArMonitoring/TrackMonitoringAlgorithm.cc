/**
 *  @file   larpandoracontent/LArMonitoring/TrackMonitoringAlgorithm.cc
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

#include "larpandoracontent/LArMonitoring/TrackMonitoringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackMonitoringAlgorithm::TrackMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackMonitoringAlgorithm::Run()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap);

    this->MatchPfoToMCParticle(PROTON, *pPfoList, mcToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackMonitoringAlgorithm::MatchPfoToMCParticle(const int pdg, const PfoList &pfoList, const LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, -1, 1));

    typedef std::tuple<const Pfo *, const CaloHitList, float> PfoIntersectionTuple;
    std::unordered_map<const MCParticle *, std::vector<PfoIntersectionTuple>> mcToPfosMap;
    for (const auto &[pMC, caloHits] : mcToHitsMap)
    {
        if (std::abs(pMC->GetParticleId()) != std::abs(pdg))
            continue;
        float totalMCAdc{0.f};
        for (const CaloHit *const pHit : caloHits)
            totalMCAdc += pHit->GetInputEnergy();
        if (totalMCAdc <= 0.f)
            continue;

        for (const auto &pPfo : pfoList)
        {
            CaloHitList pfoCaloHits;
            LArPfoHelper::GetAllCaloHits2D(pPfo, pfoCaloHits);

            CaloHitVector intersection;
            std::set_intersection(caloHits.begin(), caloHits.end(), pfoCaloHits.begin(), pfoCaloHits.end(), std::back_inserter(intersection));

            if (!intersection.empty())
            {
                float sharedAdc{0.f};
                for (const CaloHit *const pHit : intersection)
                    sharedAdc += pHit->GetInputEnergy();

                const float sharedFraction{sharedAdc / totalMCAdc};
                CaloHitList sharedHits(intersection.begin(), intersection.end());
                mcToPfosMap[pMC].emplace_back(std::make_tuple(pPfo, sharedHits, sharedFraction));
            }
        }
    }

    for (const auto &[pMC, tuples] : mcToPfosMap)
    {
        std::cout << "MCParticle (PDG: " << pMC->GetParticleId() << ") is matched to " << tuples.size() << " PFO(s):" << std::endl;
        CaloHitList mcCaloHitsU, mcCaloHitsV, mcCaloHitsW;
        for (const CaloHit *const pHit : mcToHitsMap.at(pMC))
        {
            if (pHit->GetHitType() == TPC_VIEW_U)
                mcCaloHitsU.emplace_back(pHit);
            else if (pHit->GetHitType() == TPC_VIEW_V)
                mcCaloHitsV.emplace_back(pHit);
            else if (pHit->GetHitType() == TPC_VIEW_W)
                mcCaloHitsW.emplace_back(pHit);
        }
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mcCaloHitsU, "U mc", BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mcCaloHitsV, "V mc", BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mcCaloHitsW, "W mc", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        for (const auto &[pPfo, sharedHits, sharedFraction] : tuples)
        {
            CaloHitList pfoCaloHitsU, pfoCaloHitsV, pfoCaloHitsW;
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoCaloHitsU);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoCaloHitsV);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoCaloHitsW);

            for (const CaloHit * pHit : pfoCaloHitsU)
            {
                bool isMCHit = std::find(sharedHits.begin(), sharedHits.end(), pHit) != sharedHits.end();
                std::cout << "  PFO hit (U, " << isMCHit << "): " << pHit->GetInputEnergy() << " ADC, " <<
                    pHit->GetPositionVector().GetX() << ", " << pHit->GetPositionVector().GetZ() << std::endl;
            }

            for (const CaloHit * pHit : pfoCaloHitsV)
            {
                bool isMCHit = std::find(sharedHits.begin(), sharedHits.end(), pHit) != sharedHits.end();
                std::cout << "  PFO hit (V, " << isMCHit << "): " << pHit->GetInputEnergy() << " ADC, " <<
                    pHit->GetPositionVector().GetX() << ", " << pHit->GetPositionVector().GetZ() << std::endl;
            }

            for (const CaloHit * pHit : pfoCaloHitsW)
            {
                bool isMCHit = std::find(sharedHits.begin(), sharedHits.end(), pHit) != sharedHits.end();
                std::cout << "  PFO hit (W, " << isMCHit << "): " << pHit->GetInputEnergy() << " ADC, " <<
                    pHit->GetPositionVector().GetX() << ", " << pHit->GetPositionVector().GetZ() << std::endl;
            }

            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoCaloHitsU, "U pfo", GREEN));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoCaloHitsV, "V pfo", GREEN));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoCaloHitsW, "W pfo", GREEN));
            pfoCaloHitsU.clear(); pfoCaloHitsV.clear(); pfoCaloHitsW.clear();

            for (const CaloHit *const pHit : sharedHits)
            {
                if (pHit->GetHitType() == TPC_VIEW_U)
                    pfoCaloHitsU.emplace_back(pHit);
                else if (pHit->GetHitType() == TPC_VIEW_V)
                    pfoCaloHitsV.emplace_back(pHit);
                else if (pHit->GetHitType() == TPC_VIEW_W)
                    pfoCaloHitsW.emplace_back(pHit);
            }
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoCaloHitsU, "U int", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoCaloHitsV, "V int", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoCaloHitsW, "W int", RED));

            CaloHitList allPfoHits;
            LArPfoHelper::GetAllCaloHits2D(pPfo, allPfoHits);
            std::cout << "  PFO (" << allPfoHits.size() << ") with " << sharedHits.size() << " shared hits, corresponding to " << sharedFraction * 100.f << "%" << std::endl;
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
