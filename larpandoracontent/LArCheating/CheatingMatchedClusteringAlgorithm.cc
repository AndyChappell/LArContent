/**
 *  @file   larpandoracontent/LArMonitoring/CheatingMatchedClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingMatchedClusteringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingMatchedClusteringAlgorithm::CheatingMatchedClusteringAlgorithm() :
    m_visualize{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CheatingMatchedClusteringAlgorithm::~CheatingMatchedClusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingMatchedClusteringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    }
    // Collect the hits into their respective views
    std::map<HitType, const CaloHitList *> viewHitListMap;
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (!pCaloHitList->empty())
            viewHitListMap[pCaloHitList->front()->GetHitType()] = pCaloHitList;
    }
    // Associate the hits with there respective MC particles, sort each MC's hits by x position
    typedef std::map<const MCParticle *, CaloHitVector> MCHitMap;
    std::map<HitType, MCHitMap> viewMCHitMap;
    for (const auto &[view, pCaloHitList] : viewHitListMap)
    {
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                viewMCHitMap[view][pMC].emplace_back(pCaloHit);
            }
            catch (StatusCodeException &)
            {
            }
        }
        for (auto &[pMC, hits] : viewMCHitMap[view])
            std::sort(hits.begin(), hits.end(), LArClusterHelper::SortHitsByPositionInX);
    }
    // Loop over the clusters from the seed view and find the corresponding clusters in the other view
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_seedClusterListName, pClusterList));
    for (const Cluster *const pCluster : *pClusterList)
    {
        // Collect the hits by MC and sort each MC's hits by X position (this simplifies the matching process)
        CaloHitList seedClusterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(seedClusterHits);
        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        seedClusterHits.insert(seedClusterHits.end(), isolatedHits.begin(), isolatedHits.end());
        MCHitMap seedMCHitMap;
        for (const CaloHit *pCaloHit : seedClusterHits)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                seedMCHitMap[pMC].emplace_back(pCaloHit);
            }
            catch (StatusCodeException &)
            {
            }
        }
        for (auto &[pMC, hits] : seedMCHitMap)
            std::sort(hits.begin(), hits.end(), LArClusterHelper::SortHitsByPositionInX);
        // Identify the corresponding cluster in the other views
        for (auto &[pMC, hits] : seedMCHitMap)
        {
            std::map<HitType, CaloHitVector> spanHits;
            // Get the seed cluster x span
            float minX{std::numeric_limits<float>::max()}, maxX{std::numeric_limits<float>::lowest()};
            for (const CaloHit *const pCaloHit : hits)
            {
                const float x{pCaloHit->GetPositionVector().GetX()};
                const float lo{x - 0.5f * pCaloHit->GetCellSize1()};
                const float hi{x + 0.5f * pCaloHit->GetCellSize1()};
                if (lo < minX)
                    minX = lo;
                if (hi > maxX)
                    maxX = hi;
            }
            // Get the hits in the other views that overlap with the seed span
            for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                if (view == LArClusterHelper::GetClusterHitType(pCluster))
                    continue;
                for (const CaloHit *const pCaloHit : viewMCHitMap[view][pMC])
                {
                    const float x{pCaloHit->GetPositionVector().GetX()};
                    const float lo{x - 0.5f * pCaloHit->GetCellSize1()};
                    const float hi{x + 0.5f * pCaloHit->GetCellSize1()};
                    if ((x >= minX && x <= maxX) || (minX <= lo && lo <= maxX) || (minX <= hi && hi <= maxX))
                        spanHits[view].emplace_back(pCaloHit);
                }
            }
            CaloHitList selectedHits(hits.begin(), hits.end());
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &selectedHits, "seed", BLACK));
            for (const auto &[view, overlapHits] : spanHits)
            {
                CaloHitList seedHits(overlapHits.begin(), overlapHits.end());
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &seedHits, std::to_string(static_cast<int>(view)), RED));
            }
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }

    if (m_visualize)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingMatchedClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    if (m_caloHitListNames.size() != 2)
    {
        std::cout << "CheatingMatchedClusteringAlgorithm Error: Expected two calo hit list names, received " << m_caloHitListNames.size() << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

