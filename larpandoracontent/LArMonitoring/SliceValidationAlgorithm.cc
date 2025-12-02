/**
 *  @file   larpandoracontent/LArMonitoring/SliceValidationAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include <list>
#include <unordered_set>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArMonitoring/SliceValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SliceValidationAlgorithm::SliceValidationAlgorithm() :
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SliceValidationAlgorithm::~SliceValidationAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "SliceValidationAlgorithm: Unable to write to ROOT tree" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceValidationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap, false);

    MCLeadingMap mcToLeadingMap;
    this->CreateMCToLeadingMap(mcToHitsMap, mcToLeadingMap);

    SliceHitsMap sliceToHitsMap;
    this->CreateSliceToHitsMap(mcToHitsMap, mcToLeadingMap, sliceToHitsMap);
    
    for (const auto &[pSlice, caloHits] : sliceToHitsMap)
    {
        const auto &[tpcId, pMC] = pSlice;
        std::cout << "Slice MC Particle: PDG=" << std::abs(pMC->GetParticleId()) << " E=" << pMC->GetEnergy() << " TPC=" << tpcId << " NHits=" <<
            caloHits.size() << std::endl;
    }

    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    for (const Pfo *const pPfo : *pPfoList)
    {
        CaloHitList pfoCaloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, pfoCaloHits);
        std::cout << "PFO: PDG=" << std::abs(pPfo->GetParticleId()) << " NHits=" << pfoCaloHits.size() << std::endl;
    }

    if (m_visualize)
        this->VisualizeSlices(sliceToHitsMap, *pPfoList);
    this->ValidateSlices(sliceToHitsMap, mcToLeadingMap, *pPfoList);

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
    SliceHitsMap &sliceToHitsMap) const
{
    for (const auto &[pMC, caloHits] : mcToHitsMap)
    {
        const MCParticle *const pLeadingMC{mcToLeadingMap.at(pMC)};
        // Find the TPCs with hits for this MCParticle
        std::set<unsigned int> tpcIdSet;
        for (const CaloHit *const pCaloHit : caloHits)
        {
            const LArCaloHit *const pLArCaloHit{static_cast<const LArCaloHit *>(pCaloHit)};
            if (pLArCaloHit)
                tpcIdSet.insert(pLArCaloHit->GetLArTPCVolumeId());
        }
        // Initialise the slice entries
        for (const unsigned int tpcId : tpcIdSet)
        {
            const auto key{std::make_pair(tpcId, pLeadingMC)};
            if (sliceToHitsMap.find(key) == sliceToHitsMap.end())
                sliceToHitsMap[key] = CaloHitList();
        }
        // Populate the slice entries
        for (const CaloHit *const pCaloHit : caloHits)
        {
            const LArCaloHit *const pLArCaloHit{static_cast<const LArCaloHit *>(pCaloHit)};
            if (pLArCaloHit)
            {
                const unsigned int tpcId{pLArCaloHit->GetLArTPCVolumeId()};
                sliceToHitsMap[{tpcId, pLeadingMC}].emplace_back(pCaloHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::ValidateSlices(const SliceHitsMap &mcSlices, const MCLeadingMap &mcToLeadingMap, const PfoList &recoSlices) const
{
    TrueToRecoSliceMap trueToRecoSliceMap;
    this->MatchRecoToTrueSlices(mcSlices, mcToLeadingMap, recoSlices, trueToRecoSliceMap);
    this->PopulateRootTree(trueToRecoSliceMap, mcSlices);

    if (m_visualize)
        this->VisualizeSliceMatches(trueToRecoSliceMap, mcSlices);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::MatchRecoToTrueSlices(const SliceHitsMap &mcSlices, const MCLeadingMap &mcToLeadingMap, const PfoList &recoSlices,
    TrueToRecoSliceMap &trueToRecoSliceMap) const
{
    // Loop over the reco slices and find the best matching true slice based on shared calo hits
    for (const Pfo *const pPfo : recoSlices)
    {
        CaloHitList pfoCaloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, pfoCaloHits);

        // Count the number of shared hits between this reco slice and each true slice
        std::map<const SliceIdentifier, int> mcHitCountMap;
        for (const CaloHit *const pCaloHit : pfoCaloHits)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                const MCParticle *const pLeading{mcToLeadingMap.at(pMC)};
                const LArCaloHit *const pLArCaloHit{static_cast<const LArCaloHit *>(pCaloHit)};
                if (!pLArCaloHit)
                    continue;
                const auto key{std::make_pair(pLArCaloHit->GetLArTPCVolumeId(), pLeading)};
                if (mcSlices.find(key) != mcSlices.end())
                    mcHitCountMap[key]++;
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }

        // Find the true slice with the maximum number of shared hits
        SliceIdentifier bestSlice{0, nullptr};
        int maxHits{0};
        for (const auto &[slice, hitCount] : mcHitCountMap)
        {
            if (hitCount > maxHits)
            {
                maxHits = hitCount;
                bestSlice = slice;
            }
        }

        // Add this reco slice to the map for the best matching true slice
        if (bestSlice.second)
        {
            if (trueToRecoSliceMap.find(bestSlice) == trueToRecoSliceMap.end())
                trueToRecoSliceMap[bestSlice] = PfoList();
            trueToRecoSliceMap[bestSlice].emplace_back(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::PopulateRootTree(const TrueToRecoSliceMap &trueToRecoSliceMap, const SliceHitsMap &mcSlices) const
{
    //ContingencyTable<const SliceIdentifier, const Pfo *const> cTable;
    for (const auto &[trueSlice, recoSliceList] : trueToRecoSliceMap)
    {
        const MCParticle *pMC{trueSlice.second};
        const int pdg{std::abs(pMC->GetParticleId())};
        const int isNeutrino{(pdg == NU_E || pdg == NU_MU || pdg == NU_TAU) ? 1 : 0};

        const CaloHitList &trueCaloHits{mcSlices.at(trueSlice)};
        const int trueNHits{static_cast<int>(trueCaloHits.size())};

        const int recoNSlices{static_cast<int>(recoSliceList.size())};

        for (const Pfo *const recoSlice : recoSliceList)
        {
            //cTable[trueSlice][recoSlice]++;
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "is_true_neutrino", isNeutrino));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "true_slice_pdg", pdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "true_n_hits", trueNHits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "reco_n_slices", recoNSlices));

            CaloHitList recoCaloHits;
            LArPfoHelper::GetAllCaloHits(recoSlice, recoCaloHits);
            const int recoNHits{static_cast<int>(recoCaloHits.size())};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "reco_n_hits", recoNHits));

            CaloHitList intersection;
            this->FindSetIntersection(trueCaloHits, recoCaloHits, intersection);
            const float purity{(recoNHits > 0) ? static_cast<float>(intersection.size()) / static_cast<float>(recoNHits) : 0.f};
            const float completeness{(trueNHits > 0) ? static_cast<float>(intersection.size()) / static_cast<float>(trueNHits) : 0.f};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "reco_slice_purity", purity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "reco_slice_completeness", completeness));

            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
        }
    }
    //const float ari{LArMonitoringHelper::CalcRandIndex(cTable)};
    //std::cout << "Event ARI: " << ari << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::VisualizeSlices(const SliceHitsMap &mcSlices, const PfoList &recoSlices) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));

    for (const auto &[slice, caloHitList] : mcSlices)
    {
        const MCParticle *const pMC{slice.second};
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

void SliceValidationAlgorithm::VisualizeSliceMatches(const TrueToRecoSliceMap &trueToRecoSliceMap, const SliceHitsMap &mcSlices) const
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
        const MCParticle *const pMC{trueSlice.second};
        const int pdg{std::abs(pMC->GetParticleId())};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trueCaloHitsU, "MC(U): " + std::to_string(pdg), BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trueCaloHitsV, "MC(V): " + std::to_string(pdg), BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trueCaloHitsW, "MC(W): " + std::to_string(pdg), BLUE));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::FindSetIntersection(const CaloHitList &caloHitListA, const CaloHitList &caloHitListB, CaloHitList &intersection) const
{
    std::unordered_set<const CaloHit *> setB(caloHitListB.begin(), caloHitListB.end());
    for (const CaloHit *const pCaloHit : caloHitListA)
	{
        if (setB.count(pCaloHit))
            intersection.emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
