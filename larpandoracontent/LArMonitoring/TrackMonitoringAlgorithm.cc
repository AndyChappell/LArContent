/**
 *  @file   larpandoracontent/LArMonitoring/TrackMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TrackMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArObjects/LArEventTopology.h"

using namespace pandora;

namespace lar_content
{

TrackMonitoringAlgorithm::TrackMonitoringAlgorithm() :
    m_filename("track_monitoring.root"),
    m_treename("tracks"),
    m_pfoListName("RecreatedPfos"),
    m_caloHitListName("CaloHitList2D"),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackMonitoringAlgorithm::~TrackMonitoringAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename, m_filename, "UPDATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "TrackMonitoringAlgorithm: Unable to write to ROOT tree \'" << m_treename << "\' into \'" << m_filename << "\'" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackMonitoringAlgorithm::Run()
{
    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap);
    LArMCParticleHelper::MCToPfosMap mcToPfosMap;
    LArMCParticleHelper::GetMCToPfosMap(*pPfoList, mcToPfosMap);

    this->PopulateMetrics(mcToHitsMap, mcToPfosMap);

    if (m_visualize)
        this->Visualize(mcToHitsMap, *pPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackMonitoringAlgorithm::PopulateMetrics(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const LArMCParticleHelper::MCToPfosMap &mcToPfosMap) const
{
    for (const auto &[pMC, hits] : mcToHitsMap)
    {
        if (hits.empty())
            continue;

        // Check how many of the hits are missing from any PFO in the event
        int nTotalRecoHits(0);
        for (const auto &[_, pfoList] : mcToPfosMap)
        {
            for (const Pfo *pPfo : pfoList)
            {
                int nSharedHits(LArPfoHelper::GetNSharedHits(pPfo, pMC));
                if (nSharedHits > 0)
                    nTotalRecoHits += nSharedHits;
            }
        }

        const int pdg{std::abs(pMC->GetParticleId())};
        if (pdg != PROTON && pdg != PI_PLUS && pdg != E_MINUS)
            continue;
        const LArMCParticle *const pLArMC(dynamic_cast<const LArMCParticle *>(pMC));
        if (!pLArMC)
            continue;
        if (pdg == E_MINUS)
        {
            const int process(pLArMC->GetProcess());
            // We want Michel electrons - GENIE has a bug where mu- records Michel decay as capture at rest
            // Perform more detailed checks of siblings to flag a mis-recorded Michel
            if (process != MC_PROC_DECAY && process != MC_PROC_MU_MINUS_CAPTURE_AT_REST)
                continue;
            const MCParticleList &parentList(pLArMC->GetParentList());
            if (!parentList.empty())
            {
                const MCParticle *const pParentMC(parentList.front());
                if (pParentMC)
                {
                    const int parentPdg(std::abs(pParentMC->GetParticleId()));
                    if (parentPdg != MU_MINUS && parentPdg != PI_PLUS)
                        continue;
                    const MCParticleList &siblings(pParentMC->GetDaughterList());
                    bool foundNuE(false);
                    for (const MCParticle *pSiblingMC : siblings)
                    {
                        if (pSiblingMC == pLArMC)
                            continue;
                        const int siblingPdg(std::abs(pSiblingMC->GetParticleId()));
                        if (siblingPdg == NU_E)
                        {
                            foundNuE = true;
                            break;
                        }
                    }
                    if (!foundNuE)
                        continue;
                }
            }
        }

        const CartesianPointVector &trajPoints(pLArMC->GetTrajPoints());
        if (trajPoints.size() >= 2)
        {
            float trueLength{0.f};
            for (unsigned int i = 1; i < trajPoints.size(); ++i)
            {
                const CartesianVector displacement(trajPoints[i] - trajPoints[i - 1]);
                trueLength += displacement.GetMagnitude();
            }
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "true_length", trueLength));
        }
        else
        {
            const float trueLength(0.f);
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "true_length", trueLength));
        }
        const float energy(pLArMC->GetEnergy());
        const float visibleEnergy(pLArMC->GetVisibleEnergy());
        const float momentum(pLArMC->GetMomentum().GetMagnitude());
        const bool isPrimary(LArMCParticleHelper::IsPrimary(pLArMC));
        const int nRecoParticles(mcToPfosMap.count(pLArMC) ? mcToPfosMap.at(pLArMC).size() : 0);
        int lostToParent(0), lostToChild(0), lostToSibling(0), lostToTrack(0), lostToShower(0), lost(0);

        // If an MC particle is not reconstructed, find out where the hits went (sibling, parent, child)
        if (nRecoParticles == 0)
        {
            const MCParticle *pMainLoss{nullptr}, *pMainParent{nullptr}, *pMainChild{nullptr}, *pMainSibling{nullptr};
            int nMainHitsAll{0}, nMainHits{0}, nLostHits(static_cast<int>(hits.size()));

            for (const auto &[pCurrentMC, pfoList] : mcToPfosMap)
            {
                if (pCurrentMC == pLArMC)
                    continue;
                int nSharedHitsRunning(0);
                for (const Pfo *pPfo : pfoList)
                {
                    int nSharedHits(LArPfoHelper::GetNSharedHits(pPfo, pLArMC));
                    nLostHits -= nSharedHits;
                    nSharedHitsRunning += nSharedHits;
                    if (nSharedHitsRunning > nMainHitsAll)
                    {
                        nMainHitsAll = nSharedHitsRunning;
                        pMainLoss = pCurrentMC;
                    }
                }
            }

            const MCParticleList &parentList(pLArMC->GetParentList());
            if (!parentList.empty())
            {
                const MCParticle *const pParentMC(parentList.front());
                if (pParentMC && mcToPfosMap.count(pParentMC))
                {
                    int nSharedHitsRunning(0);
                    for (const Pfo *pPfo : mcToPfosMap.at(pParentMC))
                    {
                        int nSharedHits(LArPfoHelper::GetNSharedHits(pPfo, pLArMC));
                        nSharedHitsRunning += nSharedHits;
                        if (nSharedHitsRunning > nMainHits)
                        {
                            nMainHits = nSharedHitsRunning;
                            pMainParent = pParentMC;
                        }
                    }

                    const MCParticleList &siblingList(pParentMC->GetDaughterList());
                    for (const MCParticle *pSiblingMC : siblingList)
                    {
                        if (pSiblingMC == pLArMC || !mcToPfosMap.count(pSiblingMC))
                            continue;
                        int nSharedHitsRunning(0);
                        for (const Pfo *pPfo : mcToPfosMap.at(pSiblingMC))
                        {
                            int nSharedHits(LArPfoHelper::GetNSharedHits(pPfo, pLArMC));
                            nSharedHitsRunning += nSharedHits;
                            if (nSharedHitsRunning > nMainHits)
                            {
                                nMainHits = nSharedHitsRunning;
                                pMainParent = nullptr;
                                pMainSibling = pSiblingMC;
                            }
                        }
                    }
                }
            }
            const MCParticleList &childList(pLArMC->GetDaughterList());
            for (const MCParticle *pChildMC : childList)
            {
                if (!mcToPfosMap.count(pChildMC))
                    continue;
                int nSharedHitsRunning(0);
                for (const Pfo *pPfo : mcToPfosMap.at(pChildMC))
                {
                    int nSharedHits(LArPfoHelper::GetNSharedHits(pPfo, pLArMC));
                    nSharedHitsRunning += nSharedHits;
                    if (nSharedHitsRunning > nMainHits)
                    {
                        nMainHits = nSharedHitsRunning;
                        pMainParent = nullptr;
                        pMainSibling = nullptr;
                        pMainChild = pChildMC;
                    }
                }
            }
            if (nLostHits > nMainHits && nLostHits > nMainHitsAll)
                lost = 1;
            if (nMainHitsAll > 2 && (nMainHitsAll / static_cast<float>(hits.size())) > 0.1f)
            {
                // Regardless of whether or not most hits are simply lost, record losses to other particles where the contribution is more than two
                // hits and more than 10% of the total hits
                if (nMainHitsAll > nMainHits)
                {
                    int lossPdg(std::abs(pMainLoss->GetParticleId()));
                    lostToTrack = lossPdg == PHOTON || lossPdg == E_MINUS ? 0 : 1;
                    lostToShower = lossPdg == PHOTON || lossPdg == E_MINUS ? 1 : 0;
                    pMainParent = nullptr;
                    pMainChild = nullptr;
                    pMainSibling = nullptr;
                }
                else
                {
                    if (pMainParent)
                        lostToParent = 1;
                    else if (pMainChild)
                        lostToChild = 1;
                    else if (pMainSibling)
                        lostToSibling = 1;
                }
            }
        }

        // If an MC particle has matched PFOs, compute the aggregate purity and completeness of the PFOs with respect to the MC particle
        float purity(-1.f);
        float completeness(-1.f);
        int nSharedHits(0);
        int nPfoHits(0);
        const int nMCHits(hits.size());
        if (nRecoParticles > 0)
        {
            for (const auto *const pPfo : mcToPfosMap.at(pLArMC))
            {
                CaloHitList caloHitList;
                LArPfoHelper::GetAllCaloHits2D(pPfo, caloHitList);
                nSharedHits += LArPfoHelper::GetNSharedHits(pPfo, pLArMC);
                nPfoHits += caloHitList.size();
            }
            purity = nPfoHits > 0 ? nSharedHits / static_cast<float>(nPfoHits) : 0.f;
            completeness = nSharedHits / static_cast<float>(nMCHits);
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "true_energy", energy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "true_visible_energy", visibleEnergy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "true_momentum", momentum));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "true_nhits", static_cast<int>(hits.size())));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "dropped_nhits", static_cast<int>(hits.size()) - nTotalRecoHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "pdg", pdg));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "is_primary", isPrimary ? 1 : 0));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "n_reco_particles", nRecoParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "purity", purity));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "completeness", completeness));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "lost", lost));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "lost_to_parent", lostToParent));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "lost_to_child", lostToChild));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "lost_to_sibling", lostToSibling));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "lost_to_track", lostToTrack));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "lost_to_shower", lostToShower));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackMonitoringAlgorithm::Visualize(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const pandora::PfoList &pfoList) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        CaloHitList pfoHits;
        for (const auto *const pPfo : pfoList)
            LArPfoHelper::GetCaloHits(pPfo, view, pfoHits);
        CaloHitSet pfoHitSet(pfoHits.begin(), pfoHits.end());

        for (const auto &[pMC, hits] : mcToHitsMap)
        {
            if (hits.empty())
                continue;

            const int pdg{std::abs(pMC->GetParticleId())};
            CaloHitList matchedCaloHits, droppedCaloHits;
            for (const CaloHit *const pCaloHit : hits)
            {
                if (pCaloHit->GetHitType() != view)
                    continue;
                if (pfoHitSet.find(pCaloHit) != pfoHitSet.end())
                    matchedCaloHits.emplace_back(pCaloHit);
                else
                    droppedCaloHits.emplace_back(pCaloHit);
            }
            if (!matchedCaloHits.empty())
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &matchedCaloHits, std::to_string(pdg) + " - matched", BLACK));
            }
            if (!droppedCaloHits.empty())
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &droppedCaloHits, std::to_string(pdg) + " - dropped", RED));
            }
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
