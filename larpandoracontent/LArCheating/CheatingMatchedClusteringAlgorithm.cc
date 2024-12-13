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
    std::set<const CaloHit *> usedHits;
    for (const Cluster *const pCluster : *pClusterList)
    {
        const HitType seedView{LArClusterHelper::GetClusterHitType(pCluster)};
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
        std::map<HitType, CaloHitList> matches;
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
                if (view == seedView)
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
            std::vector<HitType> views;
            std::transform(viewMCHitMap.begin(), viewMCHitMap.end(), std::back_inserter(views), [](const auto &pair){ return pair.first; });
            if (views.size() == 2)
            {
                // If we have three views populated, find the overlapping hits
                for (const CaloHit *const pCaloHit : hits)
                {
                    const float x0{pCaloHit->GetPositionVector().GetX()};
                    const float lo0{x0 - 0.5f * pCaloHit->GetCellSize1()};
                    const float hi0{x0 + 0.5f * pCaloHit->GetCellSize1()};
                    CaloHitVector view1Hits;
                    for (const CaloHit *const pCaloHit1 : viewMCHitMap[views.at(0)][pMC])
                    {
                        const float x1{pCaloHit1->GetPositionVector().GetX()};
                        const float lo1{x1 - 0.5f * pCaloHit1->GetCellSize1()};
                        const float hi1{x1 + 0.5f * pCaloHit1->GetCellSize1()};
                        if ((x0 >= lo1 && x0 <= hi1) || (x1 >= lo0 && x1 <= hi0))
                            view1Hits.emplace_back(pCaloHit1);
                    }
                    CaloHitVector view2Hits;
                    for (const CaloHit *const pCaloHit2 : viewMCHitMap[views.at(1)][pMC])
                    {
                        const float x2{pCaloHit2->GetPositionVector().GetX()};
                        const float lo2{x2 - 0.5f * pCaloHit2->GetCellSize1()};
                        const float hi2{x2 + 0.5f * pCaloHit2->GetCellSize1()};
                        if ((x0 >= lo2 && x0 <= hi2) || (x2 >= lo0 && x2 <= hi0))
                            view2Hits.emplace_back(pCaloHit2);
                    }
                    // Identify the best match, if it exists, within those candidates
                    float bestChi2{std::numeric_limits<float>::max()};
                    std::pair<const CaloHit *, const CaloHit *> match{nullptr, nullptr};
                    for (const CaloHit *const pCaloHit1 : view1Hits)
                    {
                        if (usedHits.find(pCaloHit1) != usedHits.end())
                            continue;
                        for (const CaloHit *const pCaloHit2 : view2Hits)
                        {
                            if (usedHits.find(pCaloHit2) != usedHits.end())
                                continue;
                            float chi2{0.f};
                            const CartesianVector &pos0(pCaloHit->GetPositionVector()), &pos1(pCaloHit1->GetPositionVector()),
                                &pos2(pCaloHit2->GetPositionVector());
                            CartesianVector pos3D(0, 0, 0);
                            LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), seedView, views[0], views[1], pos0, pos1, pos2, pos3D, chi2);
                            float err{LArGeometryHelper::GetSigmaUVW(this->GetPandora())};
                            chi2 *= err * err;
                            // Adjust the chi2 to factor in drift error
                            err += 0.5f * pCaloHit->GetCellSize1() + 0.5f * pCaloHit1->GetCellSize1() + 0.5f * pCaloHit2->GetCellSize1();
                            chi2 /= err * err;
                            if (chi2 < 1 && chi2 < bestChi2)
                            {
                                bestChi2 = chi2;
                                match = std::make_pair(pCaloHit1, pCaloHit2);
                            }
                        }
                    }
                    // Record the best match
                    if (bestChi2 < std::numeric_limits<float>::max())
                    {
                        matches[views[0]].emplace_back(match.first);
                        matches[views[1]].emplace_back(match.second);
                        usedHits.insert(pCaloHit);
                        usedHits.insert(match.first);
                        usedHits.insert(match.second);
                    }
                }
            }
            else if (views.size() == 1)
            {
            }
        }
        if (m_visualize)
        {
            std::vector<HitType> views;
            for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                if (view == seedView)
                    continue;
                views.emplace_back(view);
            }
            const ClusterList seedClusterList({pCluster});
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &seedClusterList, "seed", BLACK));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &matches[views[0]], std::to_string(views[0]), RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &matches[views[1]], std::to_string(views[1]), BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
    // Check for unmatched hits
    std::map<HitType, CaloHitList> unusedViewHitListMap;
    for (const auto &[view, pCaloHitList] : viewHitListMap)
    {
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            if (usedHits.find(pCaloHit) == usedHits.end())
                unusedViewHitListMap[view].emplace_back(pCaloHit);
        }
    }
    std::cout << "Unused hits: " << unusedViewHitListMap[TPC_VIEW_U].size() << " " << unusedViewHitListMap[TPC_VIEW_V].size() <<
        " " << unusedViewHitListMap[TPC_VIEW_W].size() << std::endl;

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

