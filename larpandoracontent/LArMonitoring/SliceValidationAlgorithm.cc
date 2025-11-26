/**
 *  @file   larpandoracontent/LArMonitoring/SliceValidationAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/SliceValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SliceValidationAlgorithm::SliceValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceValidationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    this->CreateMCToHitsMap(*pCaloHitList, mcToHitsMap);

    MCLeadingMap mcToLeadingMap;
    this->CreateMCToLeadingMap(mcToHitsMap, mcToLeadingMap);

    LArMCParticleHelper::MCContributionMap sliceToHitsMap;
    this->CreateSliceToHitsMap(mcToHitsMap, mcToLeadingMap, sliceToHitsMap);
    
    for (const auto &[pSliceMC, caloHits] : sliceToHitsMap)
    {
        std::cout << "Slice MC Particle: PDG=" << std::abs(pSliceMC->GetParticleId()) << " E=" << pSliceMC->GetEnergy()
                  << " NHits=" << caloHits.size() << std::endl;
    }

    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    for (const Pfo *const pPfo : *pPfoList)
    {
        CaloHitList pfoCaloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, pfoCaloHits);
        std::cout << "PFO: PDG=" << std::abs(pPfo->GetParticleId()) << " NHits=" << pfoCaloHits.size() << std::endl;
    }

    this->VisualizeSlices(sliceToHitsMap, *pPfoList);

    this->ValidateSlices(sliceToHitsMap, *pPfoList);

    // May want to veto slices where the leading particle is a neutron or photon, as these tend to be diffuse

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::CreateMCToHitsMap(const CaloHitList &caloHitList, LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            if (mcToHitsMap.find(pMC) == mcToHitsMap.end())
                mcToHitsMap[pMC] = CaloHitList();
            mcToHitsMap[pMC].emplace_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::CreateMCToLeadingMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCLeadingMap &mcToLeadingMap) const
{
    for (const auto &[pMC, _] : mcToHitsMap)
    {
        const MCParticle *pParent{pMC};
        while (!pParent->GetParentList().empty())
            pParent = pParent->GetParentList().front();
        mcToLeadingMap[pMC] = pParent;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::CreateSliceToHitsMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCLeadingMap &mcToLeadingMap,
    LArMCParticleHelper::MCContributionMap &sliceToHitsMap) const
{
    for (const auto &[pMC, caloHits] : mcToHitsMap)
    {
        const MCParticle *const pLeadingMC{mcToLeadingMap.at(pMC)};
        if (sliceToHitsMap.find(pLeadingMC) == sliceToHitsMap.end())
            sliceToHitsMap[pLeadingMC] = CaloHitList();
        sliceToHitsMap[pLeadingMC].insert(sliceToHitsMap.at(pLeadingMC).end(), caloHits.begin(), caloHits.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::ValidateSlices(const LArMCParticleHelper::MCContributionMap &mcSlices, const PfoList &recoSlices) const
{
    TrueToRecoSliceMap trueToRecoSliceMap;
    this->MatchRecoToTrueSlices(mcSlices, recoSlices, trueToRecoSliceMap);

    this->VisualizeSliceMatches(trueToRecoSliceMap, mcSlices);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::MatchRecoToTrueSlices(const LArMCParticleHelper::MCContributionMap &mcSlices, const PfoList &recoSlices,
    TrueToRecoSliceMap &trueToRecoSliceMap) const
{
    // Loop over the reco slices and find the best matching true slice based on shared calo hits
    for (const Pfo *const pPfo : recoSlices)
    {
        CaloHitList pfoCaloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, pfoCaloHits);

        // Count the number of shared hits between this reco slice and each true slice
        std::map<const MCParticle *, int> mcHitCountMap;
        for (const CaloHit *const pCaloHit : pfoCaloHits)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (mcSlices.find(pMC) != mcSlices.end())
                    mcHitCountMap[pMC]++;
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }

        // Find the true slice with the maximum number of shared hits
        const MCParticle *pBestMC{nullptr};
        int maxHits{0};
        for (const auto &[pMC, hitCount] : mcHitCountMap)
        {
            if (hitCount > maxHits)
            {
                maxHits = hitCount;
                pBestMC = pMC;
            }
        }

        // Add this reco slice to the map for the best matching true slice
        if (pBestMC)
        {
            if (trueToRecoSliceMap.find(pBestMC) == trueToRecoSliceMap.end())
                trueToRecoSliceMap[pBestMC] = PfoList();
            trueToRecoSliceMap[pBestMC].emplace_back(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::VisualizeSlices(const LArMCParticleHelper::MCContributionMap &mcSlices, const PfoList &recoSlices) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));

    for (const auto &[pMC, caloHitList] : mcSlices)
    {
        CaloHitList caloHitsW;
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (pCaloHit->GetHitType() == TPC_VIEW_W)
                caloHitsW.emplace_back(pCaloHit);
        }
        const int pdg{std::abs(pMC->GetParticleId())};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsW, "MC: " + std::to_string(pdg), BLUE));
    }

    for (const Pfo *const pPfo : recoSlices)
    {
        CaloHitList pfoCaloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, pfoCaloHits);
        CaloHitList caloHitsW;
        for (const CaloHit *const pCaloHit : pfoCaloHits)
        {
            if (pCaloHit->GetHitType() == TPC_VIEW_W)
                caloHitsW.emplace_back(pCaloHit);
        }
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsW, "Reco", RED));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::VisualizeSliceMatches(const TrueToRecoSliceMap &trueToRecoSliceMap, const LArMCParticleHelper::MCContributionMap &mcSlices) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));

    for (const auto &[trueSlice, recoSliceList] : trueToRecoSliceMap)
    {
        for (const Pfo *const recoSlice : recoSliceList)
        {
            CaloHitList recoCaloHits;
            LArPfoHelper::GetAllCaloHits(recoSlice, recoCaloHits);
            CaloHitList recoCaloHitsU, recoCaloHitsV, recoCaloHitsW;
            for (const CaloHit *const pCaloHit : recoCaloHits)
            {
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                        recoCaloHitsU.emplace_back(pCaloHit);
                        break;
                    case TPC_VIEW_V:
                        recoCaloHitsV.emplace_back(pCaloHit);
                        break;
                    case TPC_VIEW_W:
                        recoCaloHitsW.emplace_back(pCaloHit);
                        break;
                    default:
                        break;
                }
            }
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &recoCaloHitsU, "Reco(U)", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &recoCaloHitsV, "Reco(V)", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &recoCaloHitsW, "Reco(W)", RED));
        }
        CaloHitList trueCaloHits{mcSlices.at(trueSlice)};
        CaloHitList trueCaloHitsU, trueCaloHitsV, trueCaloHitsW;
        for (const CaloHit *const pCaloHit : trueCaloHits)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    trueCaloHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    trueCaloHitsV.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    trueCaloHitsW.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
        const int pdg{std::abs(trueSlice->GetParticleId())};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trueCaloHitsU, "MC(U): " + std::to_string(pdg), BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trueCaloHitsV, "MC(V): " + std::to_string(pdg), BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trueCaloHitsW, "MC(W): " + std::to_string(pdg), BLUE));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
