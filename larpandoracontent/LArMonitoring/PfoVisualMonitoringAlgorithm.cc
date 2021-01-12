/**
 *  @file   larpandoracontent/LArMonitoring/PfoVisualMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the PFO visualisation monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/PfoVisualMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

PfoVisualMonitoringAlgorithm::PfoVisualMonitoringAlgorithm() :
    m_foldToPrimaries(false),
    m_foldToShowers(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoVisualMonitoringAlgorithm::~PfoVisualMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoVisualMonitoringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 0;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_minHitSharingFraction = 0;
    parameters.m_foldBackHierarchy = m_foldToPrimaries;
    parameters.m_foldToLeadingShower = m_foldToShowers;

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

    PfoList finalStatePfos;
    for (const ParticleFlowObject *const pPfo : allConnectedPfos)
    {
        if (LArPfoHelper::IsFinalState(pPfo))
            finalStatePfos.push_back(pPfo);
    }

    // ATTN - Usually do this based on all MC hits, rather than target, but in this case the target thresholds are equivalent to all
    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, targetMCParticleToHitsMap, pfoToHitsMap, parameters.m_foldBackHierarchy,
        parameters.m_foldToLeadingShower);

    if (pfoToHitsMap.empty())
        return STATUS_CODE_SUCCESS;

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, LArMCParticleHelper::MCContributionMapVector({targetMCParticleToHitsMap}),
        pfoToMCHitSharingMap, mcToPfoHitSharingMap);

    if (pfoToMCHitSharingMap.empty())
        return STATUS_CODE_SUCCESS;

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    std::map<std::string, Color> colors = {{"mu", PINK}, {"e", GREEN}, {"gamma", ORANGE},
        {"kaon", GRAY}, {"pi", CYAN}, {"p", BLUE}, {"other", BLACK}};
    std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    std::map<std::string, CaloHitList> uHits, vHits, wHits;
    unsigned int pfoId{0};
    for (const auto [ pPfo, mcToHitsVector ] : pfoToMCHitSharingMap)
    {
        unsigned int nHitsSharedWithBestMC{0};
        const MCParticle *bestMC{nullptr};
        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcCaloHitListPair : mcToHitsVector)
        {
            const MCParticle *const pMC{mcCaloHitListPair.first};
            const CaloHitList &sharedHits{mcCaloHitListPair.second};

            if (sharedHits.size() > nHitsSharedWithBestMC)
            {
                nHitsSharedWithBestMC = sharedHits.size();
                bestMC = pMC;
            }
        }
        if (!bestMC)
            continue;
        std::string prefix("pfo_" + std::to_string(pfoId) + "_");
        for (const auto [ key, value ] : keys)
        {
            uHits[prefix + value] = CaloHitList();
            vHits[prefix + value] = CaloHitList();
            wHits[prefix + value] = CaloHitList();
        }
        uHits[prefix + "other"] = CaloHitList();
        vHits[prefix + "other"] = CaloHitList();
        wHits[prefix + "other"] = CaloHitList();

        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcCaloHitListPair : mcToHitsVector)
        {
            const CaloHitList &sharedHits{mcCaloHitListPair.second};
            for (const CaloHit *pCaloHit : sharedHits)
            {
                const HitType view{pCaloHit->GetHitType()};

                try
                {
                    const int pdg{std::abs(bestMC->GetParticleId())};
                    std::string key("other");
                    if (keys.find(pdg) != keys.end())
                        key = keys[pdg];
                    std::string fullKey{prefix + key};

                    if (view == HitType::TPC_VIEW_U)
                        uHits[fullKey].push_back(pCaloHit);
                    else if (view == HitType::TPC_VIEW_V)
                        vHits[fullKey].push_back(pCaloHit);
                    else
                        wHits[fullKey].push_back(pCaloHit);
                }
                catch (const StatusCodeException&)
                {
                    continue;
                }
            }
        }

        for (const auto [ key, value ] : keys)
        {
            if (uHits[prefix + value].size() > 0)
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits[prefix + value], prefix + value + "_u", colors[value]));
            if (vHits[prefix + value].size() > 0)
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits[prefix + value], prefix + value + "_v", colors[value]));
            if (wHits[prefix + value].size() > 0)
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits[prefix + value], prefix + value + "_w", colors[value]));
        }
        ++pfoId;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoVisualMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToShowers", m_foldToShowers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

