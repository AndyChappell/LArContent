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
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

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
    if (!pClusterList || pClusterList->empty())
        return STATUS_CODE_SUCCESS;
    std::set<const CaloHit *> usedHits;
    std::map<const Cluster*, std::map<HitType, CaloHitList>> clusterMatches;
    const HitType seedView{LArClusterHelper::GetClusterHitType(pClusterList->front())};
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
                clusterMatches[pCluster][seedView].emplace_back(pCaloHit);
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
                        clusterMatches[pCluster][views[0]].emplace_back(match.first);
                        clusterMatches[pCluster][views[1]].emplace_back(match.second);
                    }
                }
            }
            else if (views.size() == 1)
            {
                // If we have two views populated, find the overlapping hits
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
                    // Identify the best match, if it exists, within those candidates
                    float bestChi2{std::numeric_limits<float>::max()};
                    const CaloHit *match{nullptr};
                    for (const CaloHit *const pCaloHit1 : view1Hits)
                    {
                        if (usedHits.find(pCaloHit1) != usedHits.end())
                            continue;
                        float chi2{0.f};
                        const CartesianVector &pos0(pCaloHit->GetPositionVector()), &pos1(pCaloHit1->GetPositionVector());
                        CartesianVector pos3D(0, 0, 0);
                        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), seedView, views[0], pos0, pos1, pos3D, chi2);
                        float err{LArGeometryHelper::GetSigmaUVW(this->GetPandora())};
                        chi2 *= err * err;
                        // Adjust the chi2 to factor in drift error
                        err += 0.5f * pCaloHit->GetCellSize1() + 0.5f * pCaloHit1->GetCellSize1();
                        chi2 /= err * err;
                        if (chi2 < 1 && chi2 < bestChi2)
                        {
                            bestChi2 = chi2;
                            match = pCaloHit1;
                        }
                    }
                    // Record the best match
                    if (bestChi2 < std::numeric_limits<float>::max())
                    {
                        matches[views[0]].emplace_back(match);
                        usedHits.insert(pCaloHit);
                        usedHits.insert(match);
                        clusterMatches[pCluster][views[0]].emplace_back(match);
                    }
                }
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
            /*PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &seedClusterList, "seed", BLACK));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &matches[views[0]], std::to_string(views[0]), RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &matches[views[1]], std::to_string(views[1]), BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
            // Need to actually store/make these clusters
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

    for (const Cluster *const pCluster : *pClusterList)
    {
        std::map<HitType, CaloHitList> matches;
        CaloHitList seedClusterHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(seedClusterHits);

        // Organise the existing matched hits for this cluster into their corresponding views and MC particles
        viewMCHitMap.clear();
        for (const auto &[view, caloHitList] : clusterMatches[pCluster])
        {
            for (const CaloHit *pCaloHit : caloHitList)
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

        for (const auto &[view, unusedCaloHitList] : unusedViewHitListMap)
        {
            for (const CaloHit *pCaloHit : unusedCaloHitList)
            {
                if (usedHits.find(pCaloHit) != usedHits.end())
                    continue;
                try
                {
                    const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                    const CartesianVector &thisPos{pCaloHit->GetPositionVector()};
                    if (viewMCHitMap[view][pMC].size() < 2)
                        continue;
                    LArPcaHelper::WeightedPointVector points;
                    for (const CaloHit *const pClusterHit : viewMCHitMap[view][pMC])
                    {
                        const CartesianVector &otherPos{pClusterHit->GetPositionVector()};
                        const float weight{1.f / (LArGeometryHelper::GetWirePitch(this->GetPandora(), view) + thisPos.GetDistanceSquared(otherPos))};
                        points.emplace_back(std::make_pair(otherPos, weight));
                    }
                    CartesianVector centroid(0, 0, 0);
                    LArPcaHelper::EigenValues eigenValues(0, 0, 0);
                    LArPcaHelper::EigenVectors eigenVectors;
                    LArPcaHelper::RunPca(points, centroid, eigenValues, eigenVectors);
                    float minL{std::numeric_limits<float>::max()}, maxL{std::numeric_limits<float>::lowest()};
                    float minT{std::numeric_limits<float>::max()}, maxT{std::numeric_limits<float>::lowest()};
                    for (const CaloHit *const pClusterHit : viewMCHitMap[view][pMC])
                    {
                        const CartesianVector vec{pClusterHit->GetPositionVector() - centroid};
                        const float l{eigenVectors[0].GetDotProduct(vec)};
                        const float t{eigenVectors[1].GetDotProduct(vec)};
                        if (l < minL)
                            minL = l;
                        if (l > maxL)
                            maxL = l;
                        if (t < minT)
                            minT = t;
                        if (t > maxT)
                            maxT = t;
                    }
                    const CartesianVector vec{thisPos - centroid};
                    const float l{eigenVectors[0].GetDotProduct(vec)};
                    const float t{eigenVectors[1].GetDotProduct(vec)};
                    if (l >= minL && l <= maxL && t >= minT && t <= maxT)
                    {
                        matches[view].emplace_back(pCaloHit);
                        usedHits.insert(pCaloHit);
                        clusterMatches[pCluster][view].emplace_back(pCaloHit);
                    }
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
    }

    unusedViewHitListMap.clear();
    for (const auto &[view, pCaloHitList] : viewHitListMap)
    {
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            if (usedHits.find(pCaloHit) == usedHits.end())
                unusedViewHitListMap[view].emplace_back(pCaloHit);
        }
    }

    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        if (view == seedView)
            continue;
        const ClusterList *pNewClusterList{nullptr};
        std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewClusterList, temporaryListName));
        for (const auto &[pCluster, viewHitMap] : clusterMatches)
        {
            if (viewHitMap.find(view) == viewHitMap.end())
                continue;
            const Cluster *pNewCluster{nullptr};
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.insert(parameters.m_caloHitList.end(), viewHitMap.at(view).begin(), viewHitMap.at(view).end());
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
        }
        std::string clusterListName{view == TPC_VIEW_U ? "ClustersU" : view == TPC_VIEW_V ? "ClustersV" : "ClustersW"};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, clusterListName));
    }

    if (m_visualize)
    {
        for (const auto &[pCluster, viewHitMap] : clusterMatches)
        {
            for (const auto &[view, caloHitList] : viewHitMap)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "cluster_" + std::to_string(view), AUTOITER));
            }
        }
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

