/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/ConeAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the proximity association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ConeAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

bool ConeAssociationAlgorithm::Cone::IsClusterContained(const Cluster *const pCluster, const Algorithm * /*alg*/) const
{
    CartesianVector ab(m_b); ab -= m_a;
    CartesianVector ac(m_c); ac -= m_a;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    int containedHits{0}, uncontainedHits{0};
    //CaloHitList caloContained, caloUncontained;
    for (const auto layer : orderedCaloHitList)
    {
        const CaloHitList &caloHits{*layer.second};
        for (const CaloHit *pCaloHit : caloHits)
        {
            CartesianVector ap(pCaloHit->GetPositionVector()); ap -= m_a;
            if (this->IsPointContained(ab, ac, ap))
            {
                //caloContained.emplace_back(pCaloHit);
                ++containedHits;
            }
            else
            {
                //caloUncontained.emplace_back(pCaloHit);
                ++uncontainedHits;
            }
        }
    }
    /*if (containedHits)
    {
        PandoraMonitoringApi::AddLineToVisualization(alg->GetPandora(), &m_a, &m_b, "ConeA", BLUE, 1, 1);
        PandoraMonitoringApi::AddLineToVisualization(alg->GetPandora(), &m_a, &m_c, "ConeB", BLUE, 1, 1);
        PandoraMonitoringApi::VisualizeCaloHits(alg->GetPandora(), &caloContained, "In", RED);
        PandoraMonitoringApi::VisualizeCaloHits(alg->GetPandora(), &caloUncontained, "Out", BLACK);
        PandoraMonitoringApi::Pause(alg->GetPandora());
    }*/

    return containedHits && containedHits >= uncontainedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeAssociationAlgorithm::Cone::IsPointContained(const CartesianVector &ab, const CartesianVector &ac, const CartesianVector &ap) const
{
    const float dotABAB{ab.GetDotProduct(ab)};
    const float dotABAP{ab.GetDotProduct(ap)};
    const float dotACAC{ac.GetDotProduct(ac)};
    const float dotACAB{ac.GetDotProduct(ab)};
    const float dotACAP{ac.GetDotProduct(ap)};

    const float denom{dotACAC * dotABAB - dotACAB * dotACAB};
    if (std::fabs(denom) < std::numeric_limits<float>::epsilon())
        return false;
    const float invDenom{1.f / denom};
    const float u{(dotABAB * dotACAP - dotACAB * dotABAP) * invDenom};
    const float v{(dotACAC * dotABAP - dotACAB * dotACAP) * invDenom};

    return (u >= 0) && (v >= 0) && ((u + v) <= 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeAssociationAlgorithm::MatchingTriplet::MatchingTriplet(const ConeAssociationAlgorithm *const pAlgorithm,
    const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW, ClusterMergeMap &clusterMergeMap) :
    m_pAlgorithm{pAlgorithm},
    m_clustersU(clusterListU.begin(), clusterListU.end()),
    m_clustersV(clusterListV.begin(), clusterListV.end()),
    m_clustersW(clusterListW.begin(), clusterListW.end()),
    m_chi2{0.f}
{
    std::sort(m_clustersU.begin(), m_clustersU.end(), LArClusterHelper::SortByMinXPosition);
    std::sort(m_clustersV.begin(), m_clustersV.end(), LArClusterHelper::SortByMinXPosition);
    std::sort(m_clustersW.begin(), m_clustersW.end(), LArClusterHelper::SortByMinXPosition);

    this->AssessMatchQuality(clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::MatchingTriplet::AssessMatchQuality(ClusterMergeMap &clusterMergeMap)
{
    CartesianVector dummy(0.f, 0.f, 0.f);
    ClusterList mergeClustersU, mergeClustersV, mergeClustersW;
    for (const Cluster *const pClusterU : m_clustersU)
    {
        const CaloHitVector &hitsU{m_pAlgorithm->m_clusterToSortedHitsMap.at(pClusterU)};
        float xMinU{0.f}, xMaxU{0.f};
        pClusterU->GetClusterSpanX(xMinU, xMaxU);
        for (const Cluster *const pClusterV : m_clustersV)
        {
            const CaloHitVector &hitsV{m_pAlgorithm->m_clusterToSortedHitsMap.at(pClusterV)};
            float xMinV{0.f}, xMaxV{0.f};
            pClusterV->GetClusterSpanX(xMinV, xMaxV);
            for (const Cluster *const pClusterW : m_clustersW)
            {
                const CaloHitVector &hitsW{m_pAlgorithm->m_clusterToSortedHitsMap.at(pClusterW)};
                float xMinW{0.f}, xMaxW{0.f};
                pClusterW->GetClusterSpanX(xMinW, xMaxW);
                const float xMin{std::max({xMinU, xMinV, xMinW})}, xMax{std::min({xMaxU, xMaxV, xMaxW})};
                if (xMin >= xMax)
                    continue;

                // Find the hits in the overlap region
                int idxMinU{-1}, idxMaxU{-1}, idxMinV{-1}, idxMaxV{-1}, idxMinW{-1}, idxMaxW{-1};
                this->FindHitsInXRange(hitsU, xMin, xMax, idxMinU, idxMaxU);
                this->FindHitsInXRange(hitsV, xMin, xMax, idxMinV, idxMaxV);
                this->FindHitsInXRange(hitsW, xMin, xMax, idxMinW, idxMaxW);
                if (idxMinU == -1 || idxMinV == -1 || idxMinW == -1 || idxMaxU == -1 || idxMaxV == -1 || idxMaxW == -1)
                    continue;
                // Step through the hits determining the best 3D hit
                float chi2Sum{0.f};
                int nhits{0};
                for (int ui = idxMinU; ui <= idxMaxU; ++ui)
                {
                    float bestChi2{std::numeric_limits<float>::max()};
                    const float u{hitsU[ui]->GetPositionVector().GetX()};
                    bool noOverlap{false};
                    int vi{idxMinV};
                    for ( ; vi <= idxMaxV; ++vi)
                    {
                        const float v{hitsV[vi]->GetPositionVector().GetX()};
                        // Require the hits be sufficiently close in x
                        if (v < (u - 0.25f))
                            continue;
                        if (v > (u + 0.25f))
                            noOverlap = true;
                        break;
                    }
                    if (noOverlap || vi > idxMaxV)
                        continue;

                    int wi{idxMinW};
                    for ( ; wi <= idxMaxW; ++wi)
                    {
                        const float v{hitsV[vi]->GetPositionVector().GetX()};
                        const float w{hitsW[wi]->GetPositionVector().GetX()};
                        // Require the hits be sufficiently close in x
                        if (w < (u - 0.25f) || w < (v - 0.25f))
                            continue;
                        if (w > (u + 0.25f) || w > (v + 0.25f))
                            noOverlap = true;
                        break;
                    }
                    if (noOverlap || wi > idxMaxW)
                        continue;
                    const int wiRef{wi};
                    while (vi <= idxMaxV)
                    {
                        const float v{hitsV[vi]->GetPositionVector().GetX()};
                        const float w{hitsW[wi]->GetPositionVector().GetX()};
                        // ATTN: The following code acts like a pair of sliding windows that are constrained against each other. w slides
                        // 'forward' (+x) relative to v, once it no longer overlaps u and v, v slides forward relative to u, and w is reset
                        // to the start of the window and repeats. This continues until the all valid window overlaps have been assessed. 
                        if (v > u + 0.25f)
                            break;
                        if (w > u + 0.25f || w > v + 0.25f)
                        {
                            wi = wiRef;
                            ++vi;
                            continue;
                        }
                        float chi2{0.f};
                        LArGeometryHelper::MergeThreePositions3D(m_pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W,
                            hitsU[ui]->GetPositionVector(), hitsV[vi]->GetPositionVector(), hitsW[wi]->GetPositionVector(), dummy, chi2);
                        if (chi2 < bestChi2)
                            bestChi2 = chi2;
                        if (wi < idxMaxW)
                        {
                            ++wi;
                        }
                        else
                        {
                            wi = wiRef;
                            ++vi;
                        }
                    }
                    if (bestChi2 < std::numeric_limits<float>::max())
                    {
                        chi2Sum += bestChi2;
                        ++nhits;
                    }
                }
                if (nhits && (chi2Sum / nhits < .1f))
                {   // Tag clusters for merging
                    mergeClustersU.emplace_back(pClusterU);
                    mergeClustersV.emplace_back(pClusterV);
                    mergeClustersW.emplace_back(pClusterW);
                }
            }
        }
    }
    for (auto iter1 = mergeClustersU.begin(); iter1 != mergeClustersU.end(); ++iter1)
    {
        const Cluster *pCluster1{*iter1};
        for (auto iter2 = std::next(iter1); iter2 != mergeClustersU.end(); ++iter2)
        {
            const Cluster *pCluster2{*iter2};
            clusterMergeMap[pCluster1].emplace_back(pCluster2);
            clusterMergeMap[pCluster2].emplace_back(pCluster1);
        }
    }

    for (auto iter1 = mergeClustersV.begin(); iter1 != mergeClustersV.end(); ++iter1)
    {
        const Cluster *pCluster1{*iter1};
        for (auto iter2 = std::next(iter1); iter2 != mergeClustersV.end(); ++iter2)
        {
            const Cluster *pCluster2{*iter2};
            clusterMergeMap[pCluster1].emplace_back(pCluster2);
            clusterMergeMap[pCluster2].emplace_back(pCluster1);
        }
    }

    for (auto iter1 = mergeClustersW.begin(); iter1 != mergeClustersW.end(); ++iter1)
    {
        const Cluster *pCluster1{*iter1};
        for (auto iter2 = std::next(iter1); iter2 != mergeClustersW.end(); ++iter2)
        {
            const Cluster *pCluster2{*iter2};
            clusterMergeMap[pCluster1].emplace_back(pCluster2);
            clusterMergeMap[pCluster2].emplace_back(pCluster1);
        }
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeAssociationAlgorithm::MatchingTriplet::FindHitsInXRange(const CaloHitVector &hits, const float xMin, const float xMax, int &minIdx, int &maxIdx) const
{
    const int size{static_cast<int>(hits.size())};
    // Perform a quick check to see if the hits are fully contained in the range
    if (xMin <= hits.front()->GetPositionVector().GetX() && xMax >= hits.back()->GetPositionVector().GetX())
    {
        minIdx = 0; maxIdx = size - 1;
        return true;
    }
    minIdx = this->BinarySearch(hits, xMin, 0, hits.size() - 1, 0.25f);
    while (minIdx > 0)
    {   // When facing duplicate values we might not have returned the minimum index with such a value. Find it.
        // ATTN: We really want a floating point identical comparison here
        if (hits.at(minIdx - 1)->GetPositionVector().GetX() == xMin)
            --minIdx;
        else
            break;
    }
    maxIdx = this->BinarySearch(hits, xMax, 0, hits.size() - 1, 0.25f);
    while (maxIdx < size - 1)
    {   // When facing duplicate values we might not have returned the maximum index with such a value. Find it.
        // ATTN: We really want a floating point identical comparison here
        if (hits.at(maxIdx + 1)->GetPositionVector().GetX() == xMax)
            ++maxIdx;
        else
            break;
    }
    // If we find one coordinate and not the other, it implies [minIdx, hits.size()) or [0, maxIdx]
    const bool found{minIdx != -1 || maxIdx != -1};
    
    if (minIdx != -1 && maxIdx == -1)
        maxIdx = size - 1;
    else if (minIdx == -1 && maxIdx != -1)
        minIdx = 0;

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int ConeAssociationAlgorithm::MatchingTriplet::BinarySearch(const CaloHitVector &hits, const float x, const int p, const int r,
    const float tolerance) const
{
    if (p <= r)
    {
        int mid{(p + r) >> 1};
        if (hits[mid]->GetPositionVector().GetX() == x)
            return mid; // Note, still need to look for duplicate x values either to the left or right depending on if looking for min or max - hence do in caller
        else if (hits[mid]->GetPositionVector().GetX() > x)
            return this->BinarySearch(hits, x, p, mid - 1, tolerance);
        else
            return this->BinarySearch(hits, x, mid + 1, r, tolerance);
    }
    else
    {   // No exact match, but current pivot or adjacent indices might be within tolerance
        const int size{static_cast<int>(hits.size())};
        const int minIdx{p - 1 >= 0 ? p - 1 : 0};
        const int maxIdx{p + 1 < size ? p + 1 : size - 1};
        int bestIdx{-1};
        float closest{std::numeric_limits<float>::max()};
        for (int i = minIdx; i <= maxIdx; ++i)
        {
            const float pivotVal{hits[i]->GetPositionVector().GetX()};
            const float dx{std::abs(x - pivotVal)};
            if (dx <= tolerance && dx < closest)
            {
                bestIdx = i;
                closest = dx;
            }
        }
        return bestIdx;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeAssociationAlgorithm::OverlapGroup::OverlapGroup(const Cluster *pCluster1, const Cluster *pCluster2, const float overlapStart,
    const float overlapFinish) :
    m_overlapStart{overlapStart},
    m_overlapFinish{overlapFinish}
{
    m_clusters.emplace_back(pCluster1);
    m_clusters.emplace_back(pCluster2);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeAssociationAlgorithm::ViewCluster::ViewCluster(const Cluster *pClusterU, const Cluster *pClusterV, const Cluster *pClusterW,
    const float chi2, const float overlapStart, const float overlapFinish) :
    m_chi2{chi2},
    m_overlapStart{overlapStart},
    m_overlapFinish{overlapFinish},
    m_overlapSize{overlapFinish - overlapStart},
    m_minHitsInOverlapRegion{std::numeric_limits<int>::max()},
    m_maxHitsInOverlapRegion{0}
{
    m_clusters.emplace_back(pClusterU);
    m_clusters.emplace_back(pClusterV);
    m_clusters.emplace_back(pClusterW);
    // Find the minimum number of hits in a cluster in the overlap region
    for (const Cluster *pCluster : m_clusters)
    {
        int nHits{0};
        const OrderedCaloHitList &orderedCaloHits(pCluster->GetOrderedCaloHitList());
        for (auto iter = orderedCaloHits.begin(); iter != orderedCaloHits.end(); ++iter)
        {
            for (auto hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter)
            {
                float x{(*hitIter)->GetPositionVector().GetX()};
                if (x >= m_overlapStart && x < m_overlapFinish)
                    ++nHits;
            }
        }
        if (nHits < m_minHitsInOverlapRegion)
            m_minHitsInOverlapRegion = nHits;
        if (nHits > m_maxHitsInOverlapRegion)
            m_maxHitsInOverlapRegion = nHits;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::ViewCluster::GetSharedViews(const ViewCluster &other, std::vector<HitType> &sharedViews) const
{
    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        if (this->HasSharedCluster(other, view))
            sharedViews.emplace_back(view);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeAssociationAlgorithm::ViewCluster::HasSharedCluster(const ViewCluster &other) const
{
    return this->HasSharedCluster(other, TPC_VIEW_U) || this->HasSharedCluster(other, TPC_VIEW_V) ||
        this->HasSharedCluster(other, TPC_VIEW_W);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeAssociationAlgorithm::ViewCluster::HasSharedCluster(const ViewCluster &other, const HitType view) const
{
    for (const Cluster *pThisCluster : m_clusters)
    {
        if (LArClusterHelper::GetClusterHitType(pThisCluster) == view)
        {
            for (const Cluster *pOtherCluster : other.m_clusters)
            {
                if (pThisCluster == pOtherCluster)
                    return true;
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeAssociationAlgorithm::ConeAssociationAlgorithm() :
    m_minClusterLayers{1},
    m_runCount{0},
    m_visualize{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ConeAssociationAlgorithm::~ConeAssociationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(); iter != pClusterList->end(); ++iter)
    {
        const Cluster *const pCluster{*iter};

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers)
            continue;

        clusterVector.emplace_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::GetListOfSeedClusters(const ClusterVector &pClusterVector, ClusterVector &seedVector,
    ClusterVector &remnantVector) const
{
    for (const Cluster *const pCluster : pClusterVector)
    {
        if (pCluster->GetNCaloHits() < 10)
            remnantVector.emplace_back(pCluster);
        else
            seedVector.emplace_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterVector seedClusters, remnantClusters;
    this->GetListOfSeedClusters(clusterVector, seedClusters, remnantClusters);
    // Separate the clusters into their respective views
    ClusterList clusterListU, clusterListV, clusterListW;
    for (const Cluster *pCluster : seedClusters)
    {
        switch (LArClusterHelper::GetClusterHitType(pCluster))
        {
            case TPC_VIEW_U:
                clusterListU.emplace_back(pCluster);
                break;
            case TPC_VIEW_V:
                clusterListV.emplace_back(pCluster);
                break;
            case TPC_VIEW_W:
                clusterListW.emplace_back(pCluster);
                break;
            default:
                break;
        }
    }

    std::cout << "Seed clusters: " << clusterListU.size() << " " << clusterListV.size() << " " << clusterListW.size() << std::endl;
    std::cout << "Remnant clusters: " << remnantClusters.size() << std::endl;

    std::list<OverlapGroup> overlapGroupVector;
    for (const Cluster *const pClusterU : clusterListU)
    {
        float xMinU{0.f}, xMaxU{0.f};
        pClusterU->GetClusterSpanX(xMinU, xMaxU);
        const float uSpan{xMaxU - xMinU};
        for (const Cluster *const pClusterV : clusterListV)
        {
            float xMinV{0.f}, xMaxV{0.f};
            pClusterV->GetClusterSpanX(xMinV, xMaxV);
            const float vSpan{xMaxV - xMinV};
            const float xMin{std::max(xMinU, xMinV)}, xMax{std::min(xMaxU, xMaxV)};
            const float span{xMax - xMin};
            if (xMin >= xMax)
                continue;
            if (uSpan < std::numeric_limits<float>::epsilon() && vSpan < std::numeric_limits<float>::epsilon())
            {   // Both isochronous, group them
                overlapGroupVector.emplace_back(OverlapGroup(pClusterU, pClusterV, xMin, xMax));
            }
            else if (uSpan > std::numeric_limits<float>::epsilon() && (span / uSpan > 0.5f) &&
                vSpan > std::numeric_limits<float>::epsilon() && (span / vSpan > 0.5f))
            {   // Overlap region is at least half the span of each cluster, group them
                overlapGroupVector.emplace_back(OverlapGroup(pClusterU, pClusterV, xMin, xMax));
            }
        }

        for (const Cluster *const pClusterW : clusterListW)
        {
            float xMinW{0.f}, xMaxW{0.f};
            pClusterW->GetClusterSpanX(xMinW, xMaxW);
            const float wSpan{xMaxW - xMinW};
            const float xMin{std::max(xMinU, xMinW)}, xMax{std::min(xMaxU, xMaxW)};
            const float span{xMax - xMin};
            if (xMin >= xMax)
                continue;
            if (uSpan < std::numeric_limits<float>::epsilon() && wSpan < std::numeric_limits<float>::epsilon())
            {   // Both isochronous, group them
                overlapGroupVector.emplace_back(OverlapGroup(pClusterU, pClusterW, xMin, xMax));
            }
            else if (uSpan > std::numeric_limits<float>::epsilon() && (span / uSpan > 0.5f) &&
                wSpan > std::numeric_limits<float>::epsilon() && (span / wSpan > 0.5f))
            {   // Overlap region is at least half the span of each cluster, group them
                overlapGroupVector.emplace_back(OverlapGroup(pClusterU, pClusterW, xMin, xMax));
            }
        }
    }

    for (const Cluster *const pClusterV : clusterListV)
    {
        float xMinV{0.f}, xMaxV{0.f};
        pClusterV->GetClusterSpanX(xMinV, xMaxV);
        const float vSpan{xMaxV - xMinV};

        for (const Cluster *const pClusterW : clusterListW)
        {
            float xMinW{0.f}, xMaxW{0.f};
            pClusterW->GetClusterSpanX(xMinW, xMaxW);
            const float wSpan{xMaxW - xMinW};
            const float xMin{std::max(xMinV, xMinW)}, xMax{std::min(xMaxV, xMaxW)};
            const float span{xMax - xMin};
            if (xMin >= xMax)
                continue;
            if (vSpan < std::numeric_limits<float>::epsilon() && wSpan < std::numeric_limits<float>::epsilon())
            {   // Both isochronous, group them
                overlapGroupVector.emplace_back(OverlapGroup(pClusterV, pClusterW, xMin, xMax));
            }
            else if (vSpan > std::numeric_limits<float>::epsilon() && (span / vSpan > 0.5f) &&
                wSpan > std::numeric_limits<float>::epsilon() && (span / wSpan > 0.5f))
            {   // Overlap region is at least half the span of each cluster, group them
                overlapGroupVector.emplace_back(OverlapGroup(pClusterV, pClusterW, xMin, xMax));
            }
        }
    }

    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    for (OverlapGroup &group : overlapGroupVector)
    {
        const Cluster *const pCluster1{group.GetCluster1()};
        const Cluster *const pCluster2{group.GetCluster2()};
        CartesianVector start1(0.f, 0.f, 0.f), finish1(0.f, 0.f, 0.f), transverse1(0.f, 0.f, 0.f);
        float length1{0.f}, ratio1{0.f};
        this->GetConeParameters(pCluster1, start1, finish1, transverse1, length1, ratio1);
        CartesianVector start2(0.f, 0.f, 0.f), finish2(0.f, 0.f, 0.f), transverse2(0.f, 0.f, 0.f);
        float length2{0.f}, ratio2{0.f};
        this->GetConeParameters(pCluster2, start2, finish2, transverse2, length2, ratio2);

        // Make sure start is low x (or low z if x is equal)
        if (start1.GetX() > finish1.GetX())
        {
            std::swap(start1, finish1);
        }
        else if (start1.GetX() == finish1.GetX())
        {   // Check z and swap if start is downstream of finish
            if (start1.GetZ() > finish1.GetZ())
                std::swap(start1, finish1);
        }
        if (start2.GetX() > finish2.GetX())
        {
            std::swap(start2, finish2);
        }
        else if (start2.GetX() == finish2.GetX())
        {
            if (start2.GetZ() > finish2.GetZ())
                std::swap(start2, finish2);
        }

        //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start1, &finish1, "U", RED, 3, 1));
        //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start2, &finish2, "V", GREEN, 3, 1));

        // Can then try the cone in both directions and see what produces the best outcome
        // Start with a longitudinal cone, then try transverse
        CartesianVector start3(0.f, 0.f, 0.f), finish3(0.f, 0.f, 0.f);
        std::vector<HitType> views({TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W});
        views.erase(std::find(views.begin(), views.end(), LArClusterHelper::GetClusterHitType(pCluster1)));
        views.erase(std::find(views.begin(), views.end(), LArClusterHelper::GetClusterHitType(pCluster2)));

        const HitType projView{views.front()};
        std::cout << "View: " << projView << std::endl;
        switch (projView)
        {
            case TPC_VIEW_U:
            {
                this->AdjustEndPoints(start1, finish1, start2, finish2);
                auto transform{PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()};
                const float z0_3{static_cast<float>(transform->VWtoU(start1.GetZ(), start2.GetZ()))},
                    z1_3{static_cast<float>(transform->VWtoU(finish1.GetZ(), finish2.GetZ()))};
                start3.SetValues(start1.GetX(), 0.f, z0_3);
                finish3.SetValues(finish1.GetX(), 0.f, z1_3);
                CartesianVector axis3(finish3); axis3 -= start3;

                // Create the cones in the three views
                CartesianVector origin1(0.f, 0.f, 0.f), vertex1A(0.f, 0.f, 0.f), vertex1B(0.f, 0.f, 0.f);
                this->MakeCone(start1, finish1, transverse1, ratio1, origin1, vertex1A, vertex1B);
                Cone cone1(origin1, vertex1A, vertex1B);
                CartesianVector origin2(0.f, 0.f, 0.f), vertex2A(0.f, 0.f, 0.f), vertex2B(0.f, 0.f, 0.f);
                this->MakeCone(start2, finish2, transverse2, ratio2, origin2, vertex2A, vertex2B);
                Cone cone2(origin2, vertex2A, vertex2B);
                CartesianVector origin3(0.f, 0.f, 0.f), vertex3A(0.f, 0.f, 0.f), vertex3B(0.f, 0.f, 0.f);
                CartesianVector transverse3(-axis3.GetZ(), 0.f, axis3.GetX());
                transverse3 = transverse3.GetUnitVector();
                this->MakeCone(start3, finish3, transverse3, ratio1 < ratio2 ? ratio1 : ratio2, origin3, vertex3A, vertex3B);
                Cone cone3(origin3, vertex3A, vertex3B);

                this->UpdateClusterMergeMap(cone3, cone1, cone2, clusterMergeMap);

                break;
            }
            case TPC_VIEW_V:
            {
                this->AdjustEndPoints(start1, finish1, start2, finish2);
                auto transform{PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()};
                const float z0_3{static_cast<float>(transform->WUtoV(start2.GetZ(), start1.GetZ()))},
                    z1_3{static_cast<float>(transform->WUtoV(finish2.GetZ(), finish1.GetZ()))};
                start3.SetValues(start1.GetX(), 0.f, z0_3);
                finish3.SetValues(finish1.GetX(), 0.f, z1_3);
                CartesianVector axis3(finish3); axis3 -= start3;

                // Create the cones in the three views
                CartesianVector origin1(0.f, 0.f, 0.f), vertex1A(0.f, 0.f, 0.f), vertex1B(0.f, 0.f, 0.f);
                this->MakeCone(start1, finish1, transverse1, ratio1, origin1, vertex1A, vertex1B);
                Cone cone1(origin1, vertex1A, vertex1B);
                CartesianVector origin2(0.f, 0.f, 0.f), vertex2A(0.f, 0.f, 0.f), vertex2B(0.f, 0.f, 0.f);
                this->MakeCone(start2, finish2, transverse2, ratio2, origin2, vertex2A, vertex2B);
                Cone cone2(origin2, vertex2A, vertex2B);
                CartesianVector origin3(0.f, 0.f, 0.f), vertex3A(0.f, 0.f, 0.f), vertex3B(0.f, 0.f, 0.f);
                CartesianVector transverse3(-axis3.GetZ(), 0.f, axis3.GetX());
                transverse3 = transverse3.GetUnitVector();
                this->MakeCone(start3, finish3, transverse3, ratio1 < ratio2 ? ratio1 : ratio2, origin3, vertex3A, vertex3B);
                Cone cone3(origin3, vertex3A, vertex3B);

                this->UpdateClusterMergeMap(cone1, cone3, cone2, clusterMergeMap);

                break;
            }
            case TPC_VIEW_W:
            {
                this->AdjustEndPoints(start1, finish1, start2, finish2);
                auto transform{PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()};
                const float z0_3{static_cast<float>(transform->UVtoW(start1.GetZ(), start2.GetZ()))},
                    z1_3{static_cast<float>(transform->UVtoW(finish1.GetZ(), finish2.GetZ()))};
                start3.SetValues(start1.GetX(), 0.f, z0_3);
                finish3.SetValues(finish1.GetX(), 0.f, z1_3);
                CartesianVector axis3(finish3); axis3 -= start3;

                // Create the cones in the three views
                CartesianVector origin1(0.f, 0.f, 0.f), vertex1A(0.f, 0.f, 0.f), vertex1B(0.f, 0.f, 0.f);
                this->MakeCone(start1, finish1, transverse1, ratio1, origin1, vertex1A, vertex1B);
                Cone cone1(origin1, vertex1A, vertex1B);
                CartesianVector origin2(0.f, 0.f, 0.f), vertex2A(0.f, 0.f, 0.f), vertex2B(0.f, 0.f, 0.f);
                this->MakeCone(start2, finish2, transverse2, ratio2, origin2, vertex2A, vertex2B);
                Cone cone2(origin2, vertex2A, vertex2B);
                CartesianVector origin3(0.f, 0.f, 0.f), vertex3A(0.f, 0.f, 0.f), vertex3B(0.f, 0.f, 0.f);
                CartesianVector transverse3(-axis3.GetZ(), 0.f, axis3.GetX());
                transverse3 = transverse3.GetUnitVector();
                this->MakeCone(start3, finish3, transverse3, ratio1 < ratio2 ? ratio1 : ratio2, origin3, vertex3A, vertex3B);
                Cone cone3(origin3, vertex3A, vertex3B);

                this->UpdateClusterMergeMap(cone1, cone2, cone3, clusterMergeMap);

                //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start3, &finish3, "W", BLUE, 3, 1));
                //PANDORA_MONITORING_API(Pause(this->GetPandora()));
                break;
            }
            default:
                break;
        }
    }

    switch (m_runCount)
    {
        case 0:
        {
            break;
        }
        default:
        {
            break;
        }
    }
    ++m_runCount;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::AdjustEndPoints(CartesianVector &start1, CartesianVector &finish1, CartesianVector &start2,
    CartesianVector &finish2) const
{
    const float x0{std::max(start1.GetX(), start2.GetX())}, x1{std::min(finish1.GetX(), finish2.GetX())};
    float z0_1{start1.GetZ()}, z0_2{start2.GetZ()}, z1_1{finish1.GetZ()}, z1_2{finish2.GetZ()};
    const float dx1{finish1.GetX() - start1.GetX()};
    const float dx2{finish2.GetX() - start2.GetX()};
    float m1{0.f}, m2{0.f};
    if (dx1 != 0.f)
    {   // ATTN: Avoiding precisely the isochronous case, where we don't want to alter the end points
        m1 = (finish1.GetZ() - start1.GetZ()) / dx1;
        z0_1 += m1 * (x0 - start1.GetX());
        z1_1 += m1 * (x1 - finish1.GetX());
    }
    if (dx2 != 0.f)
    {   // ATTN: Avoiding precisely the isochronous case, where we don't want to alter the end points
        m2 = (finish2.GetZ() - start2.GetZ()) / dx2;
        z0_2 += m2 * (x0 - start2.GetX());
        z1_2 += m2 * (x1 - finish1.GetX());
    }
    start1.SetValues(x0, 0.f, z0_1); finish1.SetValues(x1, 0.f, z1_1);
    start2.SetValues(x0, 0.f, z0_2); finish2.SetValues(x1, 0.f, z1_2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::UpdateClusterMergeMap(const Cone &coneU, const Cone &coneV, const Cone &coneW,
    ClusterMergeMap &clusterMergeMap) const
{
    ClusterList containedU;
    for (const Cluster *pCluster : m_clustersU)
        if (coneU.IsClusterContained(pCluster, this))
            containedU.emplace_back(pCluster);

    ClusterList containedV;
    for (const Cluster *pCluster : m_clustersV)
        if (coneV.IsClusterContained(pCluster, this))
            containedV.emplace_back(pCluster);

    ClusterList containedW;
    for (const Cluster *pCluster : m_clustersW)
        if (coneW.IsClusterContained(pCluster, this))
            containedW.emplace_back(pCluster);

    MatchingTriplet triplet(this, containedU, containedV, containedW, clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::AssociateClusters(const ClusterList &clusterListU, const ClusterList &clusterListV,
const ClusterList &clusterListW, ViewClusterVector &viewClusterVector) const
{
ClusterExtremumMap clusterXminMap, clusterXmaxMap;
ClusterHitsMap clusterHitsMap;
for (const ClusterList &clusterList : {clusterListU, clusterListV, clusterListW})
        for (const Cluster *pCluster : clusterList)
            this->PopulateClusterMaps(pCluster, clusterHitsMap, clusterXminMap, clusterXmaxMap);

    for (const Cluster *pClusterW : clusterListW)
    {
        const CaloHitVector &hitsW{clusterHitsMap[pClusterW]};
        for (const Cluster *pClusterV : clusterListV)
        {
            const CaloHitVector &hitsV{clusterHitsMap[pClusterV]};
            for (const Cluster *pClusterU : clusterListU)
            {
                const CaloHitVector &hitsU{clusterHitsMap[pClusterU]};
                const float overlapStart{std::max({clusterXminMap[pClusterU], clusterXminMap[pClusterV], clusterXminMap[pClusterW]})};
                const float overlapFinish{std::min({clusterXmaxMap[pClusterU], clusterXmaxMap[pClusterV], clusterXmaxMap[pClusterW]})};
                const float overlapSize{overlapFinish - overlapStart};
                // Threshold check reduces the amount of clustering, but avoids some of the ridiculous merges
                int clustersPassingThreshold{0};
                if ((0.5f * (clusterXmaxMap[pClusterU] - clusterXminMap[pClusterU])) <= overlapSize)
                    ++clustersPassingThreshold;
                if ((0.5f * (clusterXmaxMap[pClusterV] - clusterXminMap[pClusterV])) <= overlapSize)
                    ++clustersPassingThreshold;
                if ((0.5f * (clusterXmaxMap[pClusterW] - clusterXminMap[pClusterW])) <= overlapSize)
                    ++clustersPassingThreshold;
                if (clustersPassingThreshold < 2)
                    continue;
                float meanChi2{this->AssessHitOverlap(hitsU, hitsV, hitsW, overlapStart, overlapFinish)};
                if (meanChi2 < 0.1f)
                {
                    // Then look at considering clusters together and see if PCA makes sense for each combination -
                    // need to some suitable definition of 'makes sense', especially when one cluster dominates the
                    // overlap region or dominates amongst the views, so perhaps try 3D PCA for each ViewCluster independently
                    // and see if it would extend reasonably to incorporate additional clusters from other views)
                    viewClusterVector.emplace_back(new ViewCluster(pClusterU, pClusterV, pClusterW, meanChi2, overlapStart,
                        overlapFinish));
                    ClusterList listU({pClusterU}), listV({pClusterV}), listW({pClusterW});
                    if (m_visualize)
                    {
                        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &listU, "U", RED, false));
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &listV, "V", GREEN, false));
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &listW, "W", BLUE, false));
                        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeAssociationAlgorithm::AssessHitOverlap(const pandora::CaloHitVector &hitsU, const pandora::CaloHitVector &hitsV,
    const pandora::CaloHitVector &hitsW, const float overlapStart, const float overlapFinish) const
{
    if (overlapFinish >= overlapStart)
    {   // There is an overlap, assess this triplet further
        CartesianVector dummy(0.f, 0.f, 0.f);
        float chi2Sum{0.f};
        int numMatchedBins{0};
        const int N{std::max(1, static_cast<int>(std::ceil(10.f * (overlapFinish - overlapStart))))};
        for (int bin = 0; bin < N - 1; ++bin)
        {
            CaloHitList localHitsU, localHitsV, localHitsW;
            const float start{overlapStart + 0.1f * bin}, finish{std::min(overlapStart + 0.1f * (bin + 1), overlapFinish)};
            for (const CaloHit *pCaloHit : hitsU)
            {
                const float hitX{pCaloHit->GetPositionVector().GetX()};
                if (hitX >= start && hitX < finish)
                    localHitsU.emplace_back(pCaloHit);
                if (hitX >= overlapFinish)
                    break;
            }
            for (const CaloHit *pCaloHit : hitsV)
            {
                const float hitX{pCaloHit->GetPositionVector().GetX()};
                if (hitX >= start && hitX < finish)
                    localHitsV.emplace_back(pCaloHit);
                if (hitX >= overlapFinish)
                    break;
            }
            for (const CaloHit *pCaloHit : hitsW)
            {
                const float hitX{pCaloHit->GetPositionVector().GetX()};
                if (hitX >= start && hitX < finish)
                    localHitsW.emplace_back(pCaloHit);
                if (hitX >= overlapFinish)
                    break;
            }
            float minChi2{std::numeric_limits<float>::max()};
            for (const CaloHit *pCaloHitU : localHitsU)
            {
                for (const CaloHit *pCaloHitV : localHitsV)
                {
                    for (const CaloHit *pCaloHitW : localHitsW)
                    {
                        float chi2{0.f};
                        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W,
                            pCaloHitU->GetPositionVector(), pCaloHitV->GetPositionVector(), pCaloHitW->GetPositionVector(),
                            dummy, chi2);
                        if (chi2 < minChi2)
                            minChi2 = chi2;
                    }
                }
            }
            if (minChi2 < std::numeric_limits<float>::max())
            {
                chi2Sum += minChi2;
                ++numMatchedBins;
            }
        }
        return numMatchedBins > 0 ? chi2Sum / numMatchedBins : std::numeric_limits<float>::max();
    }

    return std::numeric_limits<float>::max();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::PopulateClusterMaps(const pandora::Cluster *pCluster, ClusterHitsMap &clusterHitsMap,
    ClusterExtremumMap &clusterXminMap, ClusterExtremumMap &clusterXmaxMap) const
{
    CaloHitList hits;
    const OrderedCaloHitList &orderedCaloHits(pCluster->GetOrderedCaloHitList());
    orderedCaloHits.FillCaloHitList(hits);
    CaloHitVector sortedHits(hits.begin(), hits.end());
    std::sort(sortedHits.begin(), sortedHits.end(), LArClusterHelper::SortHitsByPositionInX);
    clusterXminMap[pCluster] = sortedHits.front()->GetPositionVector().GetX();
    clusterXmaxMap[pCluster] = sortedHits.back()->GetPositionVector().GetX();
    clusterHitsMap[pCluster] = sortedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::ConsolidateClusters(ViewClusterVector &viewClusterVector, ClusterMergeMap &clusterMergeMap) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    
    std::sort(viewClusterVector.begin(), viewClusterVector.end(), ViewCluster::Sort);
    ClusterExtremumMap clusterXminMap, clusterXmaxMap;
    ClusterHitsMap clusterHitsMap;
    for (const ViewCluster *viewCluster : viewClusterVector)
        for (const Cluster *pCluster : viewCluster->GetClusterList())
            this->PopulateClusterMaps(pCluster, clusterHitsMap, clusterXminMap, clusterXmaxMap);

    for (auto iter1 = viewClusterVector.begin(); iter1 != viewClusterVector.end(); ++iter1)
    {
        const ViewCluster viewCluster1(*(*iter1));
        const ClusterList clusters1{viewCluster1.GetClusterList()};
        for (auto iter2 = std::next(iter1); iter2 != viewClusterVector.end(); ++iter2)
        {
            const ViewCluster viewCluster2(*(*iter2));
            if (viewCluster1.HasSharedCluster(viewCluster2))
            {   // We have a shared cluster, see if we can combine additional clusters
                std::vector<HitType> sharedViews;
                viewCluster1.GetSharedViews(viewCluster2, sharedViews);
                CaloHitList hitListU, hitListV, hitListW;
                for (const pandora::Cluster *pCluster : clusters1)
                {
                    const OrderedCaloHitList &orderedCaloHits(pCluster->GetOrderedCaloHitList());

                    const HitType view(LArClusterHelper::GetClusterHitType(pCluster));
                    if (view == TPC_VIEW_U)
                        orderedCaloHits.FillCaloHitList(hitListU);
                    else if (view == TPC_VIEW_V)
                        orderedCaloHits.FillCaloHitList(hitListV);
                    else
                        orderedCaloHits.FillCaloHitList(hitListW);
                    pandora::ClusterList temp({pCluster});
                    if (m_visualize)
                    {
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &temp, std::to_string(view), BLACK, false));
                    }
                }
                const ClusterList clusters2{viewCluster2.GetClusterList()};
                for (const pandora::Cluster *pCluster : clusters2)
                {
                    if (std::find(sharedViews.begin(), sharedViews.end(), LArClusterHelper::GetClusterHitType(pCluster)) !=
                        sharedViews.end())
                        continue;
                    const OrderedCaloHitList &orderedCaloHits(pCluster->GetOrderedCaloHitList());

                    const HitType view(LArClusterHelper::GetClusterHitType(pCluster));
                    if (view == TPC_VIEW_U)
                        orderedCaloHits.FillCaloHitList(hitListU);
                    else if (view == TPC_VIEW_V)
                        orderedCaloHits.FillCaloHitList(hitListV);
                    else
                        orderedCaloHits.FillCaloHitList(hitListW);
                    pandora::ClusterList temp({pCluster});
                    if (m_visualize)
                    {
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &temp, std::to_string(view), CYAN, false));
                    }
                }
            
                CaloHitVector sortedHitsU(hitListU.begin(), hitListU.end());
                std::sort(sortedHitsU.begin(), sortedHitsU.end(), LArClusterHelper::SortHitsByPositionInX);
                CaloHitVector sortedHitsV(hitListV.begin(), hitListV.end());
                std::sort(sortedHitsV.begin(), sortedHitsV.end(), LArClusterHelper::SortHitsByPositionInX);
                CaloHitVector sortedHitsW(hitListW.begin(), hitListW.end());
                std::sort(sortedHitsW.begin(), sortedHitsW.end(), LArClusterHelper::SortHitsByPositionInX);
                const float overlapStart{std::max({
                    !sortedHitsU.empty() ? sortedHitsU.front()->GetPositionVector().GetX() : 0,
                    !sortedHitsV.empty() ? sortedHitsV.front()->GetPositionVector().GetX() : 0,
                    !sortedHitsW.empty() ? sortedHitsW.front()->GetPositionVector().GetX() : 0})};
                const float overlapFinish{std::min({
                    !sortedHitsU.empty() ? sortedHitsU.back()->GetPositionVector().GetX() : std::numeric_limits<float>::max(),
                    !sortedHitsV.empty() ? sortedHitsV.back()->GetPositionVector().GetX() : std::numeric_limits<float>::max(),
                    !sortedHitsW.empty() ? sortedHitsW.back()->GetPositionVector().GetX() : std::numeric_limits<float>::max()})};

                float chi2{this->AssessHitOverlap(sortedHitsU, sortedHitsV, sortedHitsW, overlapStart, overlapFinish)};
                const float threshold{1.1f * std::max(viewCluster1.GetChi2(), viewCluster2.GetChi2())};
                if (chi2 < threshold)
                {   // Merging these clusters would still yield a good chi2
                    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
                    {
                        if (std::find(sharedViews.begin(), sharedViews.end(), view) == sharedViews.end())
                        {
                            const Cluster *pCluster1{viewCluster1.GetCluster(view)};
                            const Cluster *pCluster2{viewCluster2.GetCluster(view)};
                            // clusterMergeMap.find(pCluster) == clusterMergeMap.end()
                            clusterMergeMap[pCluster1].emplace_back(pCluster2);
                            clusterMergeMap[pCluster2].emplace_back(pCluster1);
                        }
                    }
                }
                if (m_visualize)
                {
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitListU, "U", RED));
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitListV, "V", GREEN));
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitListW, "W", BLUE));
                    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::GrowClusters(ViewClusterVector &viewClusterVector, ClusterMergeMap &clusterMergeMap) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    
    std::sort(viewClusterVector.begin(), viewClusterVector.end(), ViewCluster::SortMax);
    ClusterExtremumMap clusterXminMap, clusterXmaxMap;
    ClusterHitsMap clusterHitsMap;
    for (const ViewCluster *viewCluster : viewClusterVector)
        for (const Cluster *pCluster : viewCluster->GetClusterList())
            this->PopulateClusterMaps(pCluster, clusterHitsMap, clusterXminMap, clusterXmaxMap);

    for (auto iter1 = viewClusterVector.begin(); iter1 != viewClusterVector.end(); ++iter1)
    {
        const ViewCluster viewCluster(*(*iter1));
        // Build a 3D PCA-based cone around this candidate set of clusters and see if it can be grown
        // Note, here it's much more important that the relative quality of the cluster association is
        // considered, as there will be conflicting cones, and we should only pick the best
        // Likely to be a multi-iteration process
        (void)clusterMergeMap;

        // Identify the two views with the most hits
        const Cluster *pClusterU{viewCluster.GetCluster(TPC_VIEW_U)};
        const Cluster *pClusterV{viewCluster.GetCluster(TPC_VIEW_V)};
        const Cluster *pClusterW{viewCluster.GetCluster(TPC_VIEW_W)};
        CaloHitList caloHitsU, caloHitsV, caloHitsW;
        OrderedCaloHitList orderedCaloHitList(pClusterU->GetOrderedCaloHitList());
        orderedCaloHitList.FillCaloHitList(caloHitsU);
        orderedCaloHitList = pClusterV->GetOrderedCaloHitList();
        orderedCaloHitList.FillCaloHitList(caloHitsV);
        orderedCaloHitList = pClusterW->GetOrderedCaloHitList();
        orderedCaloHitList.FillCaloHitList(caloHitsW);
        HitTypeVector views;
        if (caloHitsU.size() > caloHitsV.size())
        {
            views.emplace_back(TPC_VIEW_U);
            if (caloHitsV.size() > caloHitsW.size())
                views.emplace_back(TPC_VIEW_V);
            else
                views.emplace_back(TPC_VIEW_W);
        }
        else if (caloHitsU.size() > caloHitsW.size())
        {
            views.emplace_back(TPC_VIEW_U);
            views.emplace_back(TPC_VIEW_V);   
        }
        else
        {
            views.emplace_back(TPC_VIEW_V);
            views.emplace_back(TPC_VIEW_W);
        }

        // Retrieve the hits for PCA
        // Aim here is to get the PCA axes of the two 'best' views, consoldate them onto the same x range, project the end points
        // into the third view, then use the broadest transverse range to define the cone in 3D, before projecting back into the three
        // views. Check this makes sense, then see about mopping up clusters falling in this cone and some extensions of it

        const Cluster *pCluster1{viewCluster.GetCluster(views[0])};
        CartesianVector pca1Start(0.f, 0.f, 0.f), pca1Finish(0.f, 0.f, 0.f), pca1Transverse(0.f, 0.f, 0.f);
        float length1{0.f}, ratio1{1.f};
        this->GetConeParameters(pCluster1, pca1Start, pca1Finish, pca1Transverse, length1, ratio1);
        const Cluster *pCluster2{viewCluster.GetCluster(views[1])};
        CartesianVector pca2Start(0.f, 0.f, 0.f), pca2Finish(0.f, 0.f, 0.f), pca2Transverse(0.f, 0.f, 0.f);
        float length2{0.f}, ratio2{1.f};
        this->GetConeParameters(pCluster2, pca2Start, pca2Finish, pca2Transverse, length2, ratio2);
        // Get the remaining view cluster for later comparison purposes
        const Cluster *pClusterX{viewCluster.GetCluster(views[0] != TPC_VIEW_U ? TPC_VIEW_U : views[1] == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_V)};
        CartesianVector pcaXStart(0.f, 0.f, 0.f), pcaXFinish(0.f, 0.f, 0.f), pcaXTransverse(0.f, 0.f, 0.f);
        float lengthX{0.f}, ratioX{1.f};
        this->GetConeParameters(pClusterX, pcaXStart, pcaXFinish, pcaXTransverse, lengthX, ratioX);

        CartesianVector aStart(pca1Start), aFinish(pca1Finish), bStart(pca2Start), bFinish(pca2Finish);

        // Ensure PCA axes have same x starting point
        const float minX{std::min({pca1Start.GetX(), pca1Finish.GetX(), pca2Start.GetX(), pca2Finish.GetX()})};
        const float maxX{std::max({pca1Start.GetX(), pca1Finish.GetX(), pca2Start.GetX(), pca2Finish.GetX()})};

        if (pca1Start.GetX() <= pca1Finish.GetX())
        {   // PCA axis 1 moves along +x from start to finish
            const float m{pca1Finish.GetX() != pca1Start.GetX() ? (pca1Finish.GetZ() - pca1Start.GetZ()) / (pca1Finish.GetX() - pca1Start.GetX()) : 0.f};
            const float startZp{pca1Start.GetZ() + m * (minX - pca1Start.GetX())};
            pca1Start.SetValues(minX, 0.f, startZp);
            const float finishZp{pca1Finish.GetZ() + m * (maxX - pca1Finish.GetX())};
            pca1Finish.SetValues(maxX, 0.f, finishZp);
        }
        else
        {   // PCA axis 1 moves along -x from start to finish
            const float m{pca1Finish.GetX() != pca1Start.GetX() ? (pca1Finish.GetZ() - pca1Start.GetZ()) / (pca1Finish.GetX() - pca1Start.GetX()) : 0.f};
            const float startZp{pca1Start.GetZ() + m * (maxX - pca1Start.GetX())};
            pca1Start.SetValues(maxX, 0.f, startZp);
            const float finishZp{pca1Finish.GetZ() + m * (minX - pca1Finish.GetX())};
            pca1Finish.SetValues(minX, 0.f, finishZp);
        }
        if (pca2Start.GetX() <= pca2Finish.GetX())
        {   // PCA axis 2 moves along +x from start to finish
            const float m{pca2Finish.GetX() != pca2Start.GetX() ? (pca2Finish.GetZ() - pca2Start.GetZ()) / (pca2Finish.GetX() - pca2Start.GetX()) : 0.f};
            const float startZp{pca2Start.GetZ() + m * (minX - pca2Start.GetX())};
            pca2Start.SetValues(minX, 0.f, startZp);
            const float finishZp{pca2Finish.GetZ() + m * (maxX - pca2Finish.GetX())};
            pca2Finish.SetValues(maxX, 0.f, finishZp);
        }
        else
        {   // PCA axis 2 moves along -x from start to finish
            const float m{pca2Finish.GetX() != pca2Start.GetX() ? (pca2Finish.GetZ() - pca2Start.GetZ()) / (pca2Finish.GetX() - pca2Start.GetX()) : 0.f};
            const float startZp{pca2Start.GetZ() + m * (maxX - pca2Start.GetX())};
            pca2Start.SetValues(maxX, 0.f, startZp);
            const float finishZp{pca2Finish.GetZ() + m * (minX - pca2Finish.GetX())};
            pca2Finish.SetValues(minX, 0.f, finishZp);
        }

        //std::cout << "Start 1: (" << pca1Start.GetX() << "," << pca1Start.GetZ() << ") 2: (" << pca2Start.GetX() << "," << pca2Start.GetZ() <<
        //    ")" << std::endl;
        //std::cout << "Finish 1: (" << pca1Finish.GetX() << "," << pca1Finish.GetZ() << ") 2: (" << pca2Finish.GetX() << "," << pca2Finish.GetZ() <<
        //    ")" << std::endl;

        //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca1Start, &pca1Finish, "1", ORANGE, 3, 1));
        //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca2Start, &pca2Finish, "2", MAGENTA, 3, 1));

        if (pca1Start.GetX() != pca1Finish.GetX() && pca2Start.GetX() != pca2Finish.GetX())
        {   // Neither axis is isochronous, endpoints can be matched on x alone
            if (pca1Start.GetX() == pca2Start.GetX())
            {   // Matching start to start and finish to finish
                CartesianVector pca3Start(0.f, 0.f, 0.f), pca3Finish(0.f, 0.f, 0.f);
                HitType projView{TPC_VIEW_W};
                if (views[0] == TPC_VIEW_U)
                {
                    if (views[1] == TPC_VIEW_V)
                    {   // UV inputs
                        //std::cout << "UVW" << std::endl;
                        const float us{pca1Start.GetZ()}, vs{pca2Start.GetZ()};
                        const float ws{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(us, vs))};
                        pca3Start.SetValues(pca1Start.GetX(), 0.f, ws);
                        const float uf{pca1Finish.GetZ()}, vf{pca2Finish.GetZ()};
                        const float wf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(uf, vf))};
                        pca3Finish.SetValues(pca1Finish.GetX(), 0.f, wf);
                    }
                    else
                    {   // UW inputs
                        //std::cout << "WUV" << std::endl;
                        projView = TPC_VIEW_V;
                        const float us{pca1Start.GetZ()}, ws{pca2Start.GetZ()};
                        const float vs{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(ws, us))};
                        pca3Start.SetValues(pca1Start.GetX(), 0.f, vs);
                        const float uf{pca1Finish.GetZ()}, wf{pca2Finish.GetZ()};
                        const float vf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(wf, uf))};
                        pca3Finish.SetValues(pca1Finish.GetX(), 0.f, vf);
                    }
                }
                else
                {   // VW inputs
                    //std::cout << "VWU" << std::endl;
                    projView = TPC_VIEW_U;
                    const float vs{pca1Start.GetZ()}, ws{pca2Start.GetZ()};
                    const float us{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vs, ws))};
                    pca3Start.SetValues(pca1Start.GetX(), 0.f, us);
                    const float vf{pca1Finish.GetZ()}, wf{pca2Finish.GetZ()};
                    const float uf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vf, wf))};
                    pca3Finish.SetValues(pca1Finish.GetX(), 0.f, uf);
                }
                
                //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjV", BLACK, 3, 1));
                CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                //std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                if (cosOpenAngle > 0.98f)
                {
                    CartesianVector origin1(0.f, 0.f, 0.f), vertex1A(0.f, 0.f, 0.f), vertex1B(0.f, 0.f, 0.f);
                    this->MakeCone(pca1Start, pca1Finish, pca1Transverse, ratio1, origin1, vertex1A, vertex1B);
                    Cone cone1(origin1, vertex1A, vertex1B);
                    CartesianVector origin2(0.f, 0.f, 0.f), vertex2A(0.f, 0.f, 0.f), vertex2B(0.f, 0.f, 0.f);
                    this->MakeCone(pca2Start, pca2Finish, pca2Transverse, ratio2, origin2, vertex2A, vertex2B);
                    Cone cone2(origin2, vertex2A, vertex2B);
                    CartesianVector origin3(0.f, 0.f, 0.f), vertex3A(0.f, 0.f, 0.f), vertex3B(0.f, 0.f, 0.f);
                    CartesianVector projTransverse(-axis1.GetZ(), 0.f, axis1.GetX());
                    projTransverse = projTransverse.GetUnitVector();
                    if (pca3Start.GetZ() > pca3Finish.GetZ())
                        std::swap(pca3Start, pca3Finish);
                    this->MakeCone(pca3Start, pca3Finish, projTransverse, ratio1 < ratio2 ? ratio1 : ratio2, origin3, vertex3A, vertex3B);
                    Cone cone3(origin3, vertex3A, vertex3B);

                    if (projView == TPC_VIEW_W)
                    {   // Cone 1 in U, Cone 2 in V, Cone 3 in W
                        ClusterList containedU;
                        for (const Cluster *pCluster : m_clustersU)
                            if (cone1.IsClusterContained(pCluster, this))
                                containedU.emplace_back(pCluster);

                        ClusterList containedV;
                        for (const Cluster *pCluster : m_clustersV)
                            if (cone2.IsClusterContained(pCluster, this))
                                containedV.emplace_back(pCluster);

                        ClusterList containedW;
                        for (const Cluster *pCluster : m_clustersW)
                            if (cone3.IsClusterContained(pCluster, this))
                                containedW.emplace_back(pCluster);

                        //std::cout << "Contained clusters " << containedU.size() << " " << containedV.size() << " " << containedW.size() << std::endl;

                        MatchingTriplet triplet(this, containedU, containedV, containedW, clusterMergeMap);
                    }
                    else if (projView == TPC_VIEW_V)
                    {   // Cone 1 in U, Cone 2 in W, Cone 3 in V
                        ClusterList containedU;
                        for (const Cluster *pCluster : m_clustersU)
                            if (cone1.IsClusterContained(pCluster, this))
                                containedU.emplace_back(pCluster);

                        ClusterList containedV;
                        for (const Cluster *pCluster : m_clustersV)
                            if (cone3.IsClusterContained(pCluster, this))
                                containedV.emplace_back(pCluster);

                        ClusterList containedW;
                        for (const Cluster *pCluster : m_clustersW)
                            if (cone2.IsClusterContained(pCluster, this))
                                containedW.emplace_back(pCluster);

                        //std::cout << "Contained clusters " << containedU.size() << " " << containedV.size() << " " << containedW.size() << std::endl;

                        MatchingTriplet triplet(this, containedU, containedV, containedW, clusterMergeMap);
                    }
                    else
                    {   // Cone 1 in V, Cone 2 in W, Cone 3 in U
                        ClusterList containedU;
                        for (const Cluster *pCluster : m_clustersU)
                            if (cone3.IsClusterContained(pCluster, this))
                                containedU.emplace_back(pCluster);

                        ClusterList containedV;
                        for (const Cluster *pCluster : m_clustersV)
                            if (cone1.IsClusterContained(pCluster, this))
                                containedV.emplace_back(pCluster);

                        ClusterList containedW;
                        for (const Cluster *pCluster : m_clustersW)
                            if (cone2.IsClusterContained(pCluster, this))
                                containedW.emplace_back(pCluster);

                        //std::cout << "Contained clusters " << containedU.size() << " " << containedV.size() << " " << containedW.size() << std::endl;
                        MatchingTriplet triplet(this, containedU, containedV, containedW, clusterMergeMap);
                    }
                }
            }
            else
            {   // Matching start1 to finish2 and start2 to finish1
                CartesianVector pca3Start(0.f, 0.f, 0.f), pca3Finish(0.f, 0.f, 0.f);
                HitType projView{TPC_VIEW_W};
                if (views[0] == TPC_VIEW_U)
                {
                    if (views[1] == TPC_VIEW_V)
                    {   // UV inputs
                        //std::cout << "UVW" << std::endl;
                        const float us{pca1Start.GetZ()}, vs{pca2Finish.GetZ()};
                        const float ws{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(us, vs))};
                        pca3Start.SetValues(pca1Start.GetX(), 0.f, ws);
                        const float uf{pca1Finish.GetZ()}, vf{pca2Start.GetZ()};
                        const float wf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(uf, vf))};
                        pca3Finish.SetValues(pca1Finish.GetX(), 0.f, wf);
                    }
                    else
                    {   // UW inputs
                        //std::cout << "WUV" << std::endl;
                        projView = TPC_VIEW_V;
                        const float us{pca1Start.GetZ()}, ws{pca2Finish.GetZ()};
                        const float vs{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(ws, us))};
                        pca3Start.SetValues(pca1Start.GetX(), 0.f, vs);
                        const float uf{pca1Finish.GetZ()}, wf{pca2Start.GetZ()};
                        const float vf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(wf, uf))};
                        pca3Finish.SetValues(pca1Finish.GetX(), 0.f, vf);
                    }
                }
                else
                {   // VW inputs
                    //std::cout << "VWU" << std::endl;
                    projView = TPC_VIEW_U;
                    const float vs{pca1Start.GetZ()}, ws{pca2Finish.GetZ()};
                    const float us{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vs, ws))};
                    pca3Start.SetValues(pca1Start.GetX(), 0.f, us);
                    const float vf{pca1Finish.GetZ()}, wf{pca2Start.GetZ()};
                    const float uf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vf, wf))};
                    pca3Finish.SetValues(pca1Finish.GetX(), 0.f, uf);
                }

                //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjV", BLACK, 3, 1));
                CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                //std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                if (cosOpenAngle > 0.98f)
                {
                    CartesianVector origin1(0.f, 0.f, 0.f), vertex1A(0.f, 0.f, 0.f), vertex1B(0.f, 0.f, 0.f);
                    this->MakeCone(pca1Start, pca1Finish, pca1Transverse, ratio1, origin1, vertex1A, vertex1B);
                    Cone cone1(origin1, vertex1A, vertex1B);
                    CartesianVector origin2(0.f, 0.f, 0.f), vertex2A(0.f, 0.f, 0.f), vertex2B(0.f, 0.f, 0.f);
                    this->MakeCone(pca2Start, pca2Finish, pca2Transverse, ratio2, origin2, vertex2A, vertex2B);
                    Cone cone2(origin2, vertex2A, vertex2B);
                    CartesianVector origin3(0.f, 0.f, 0.f), vertex3A(0.f, 0.f, 0.f), vertex3B(0.f, 0.f, 0.f);
                    CartesianVector projTransverse(-axis1.GetZ(), 0.f, axis1.GetX());
                    projTransverse = projTransverse.GetUnitVector();
                    if (pca3Start.GetZ() > pca3Finish.GetZ())
                        std::swap(pca3Start, pca3Finish);
                    this->MakeCone(pca3Start, pca3Finish, projTransverse, ratio1 < ratio2 ? ratio1 : ratio2, origin3, vertex3A, vertex3B);
                    Cone cone3(origin3, vertex3A, vertex3B);

                    if (projView == TPC_VIEW_W)
                    {   // Cone 1 in U, Cone 2 in V, Cone 3 in W
                        ClusterList containedU;
                        for (const Cluster *pCluster : m_clustersU)
                            if (cone1.IsClusterContained(pCluster, this))
                                containedU.emplace_back(pCluster);

                        ClusterList containedV;
                        for (const Cluster *pCluster : m_clustersV)
                            if (cone2.IsClusterContained(pCluster, this))
                                containedV.emplace_back(pCluster);

                        ClusterList containedW;
                        for (const Cluster *pCluster : m_clustersW)
                            if (cone3.IsClusterContained(pCluster, this))
                                containedW.emplace_back(pCluster);

                        //std::cout << "Contained clusters " << containedU.size() << " " << containedV.size() << " " << containedW.size() << std::endl;
                        MatchingTriplet triplet(this, containedU, containedV, containedW, clusterMergeMap);
                    }
                    else if (projView == TPC_VIEW_V)
                    {   // Cone 1 in U, Cone 2 in W, Cone 3 in V
                        ClusterList containedU;
                        for (const Cluster *pCluster : m_clustersU)
                            if (cone1.IsClusterContained(pCluster, this))
                                containedU.emplace_back(pCluster);

                        ClusterList containedV;
                        for (const Cluster *pCluster : m_clustersV)
                            if (cone3.IsClusterContained(pCluster, this))
                                containedV.emplace_back(pCluster);

                        ClusterList containedW;
                        for (const Cluster *pCluster : m_clustersW)
                            if (cone2.IsClusterContained(pCluster, this))
                                containedW.emplace_back(pCluster);

                        //std::cout << "Contained clusters " << containedU.size() << " " << containedV.size() << " " << containedW.size() << std::endl;
                        MatchingTriplet triplet(this, containedU, containedV, containedW, clusterMergeMap);
                    }
                    else
                    {   // Cone 1 in V, Cone 2 in W, Cone 3 in U
                        ClusterList containedU;
                        for (const Cluster *pCluster : m_clustersU)
                            if (cone3.IsClusterContained(pCluster, this))
                                containedU.emplace_back(pCluster);

                        ClusterList containedV;
                        for (const Cluster *pCluster : m_clustersV)
                            if (cone1.IsClusterContained(pCluster, this))
                                containedV.emplace_back(pCluster);

                        ClusterList containedW;
                        for (const Cluster *pCluster : m_clustersW)
                            if (cone2.IsClusterContained(pCluster, this))
                                containedW.emplace_back(pCluster);

                        //std::cout << "Contained clusters " << containedU.size() << " " << containedV.size() << " " << containedW.size() << std::endl;
                        MatchingTriplet triplet(this, containedU, containedV, containedW, clusterMergeMap);
                    }
                }
            }
        }
        else
        {   // At least one axis is isochronous, need to check all possible endpoint matches
            std::cout << "Isochronous case" << std::endl;
        }

//        PANDORA_MONITORING_API(Pause(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::GetConeParameters(const Cluster *pCluster, CartesianVector &p0, CartesianVector &p1,
    CartesianVector &transverse, float &length, float &ratio) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);

    // Get the eigen vectors for this cluster
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(caloHits, centroid, eigenValues, eigenVectors);

    // Project the extremal hits of the cluster onto the principal axis
    LArClusterHelper::GetExtremalCoordinates(pCluster, p0, p1);
    p0 -= centroid;
    p1 -= centroid;
    const CartesianVector &principalAxis(eigenVectors[0]);
    const float vMagSquared{principalAxis.GetMagnitudeSquared()};
    const float normDot0p{p0.GetDotProduct(principalAxis) / vMagSquared};
    const float normDot1p{p1.GetDotProduct(principalAxis) / vMagSquared};
    p0.SetValues(normDot0p * principalAxis.GetX() + centroid.GetX(), normDot0p * principalAxis.GetY() + centroid.GetY(),
        normDot0p * principalAxis.GetZ() + centroid.GetZ());
    p1.SetValues(normDot1p * principalAxis.GetX() + centroid.GetX(), normDot1p * principalAxis.GetY() + centroid.GetY(),
        normDot1p * principalAxis.GetZ() + centroid.GetZ());

    // Get the length of the cluster along the principal axis
    CartesianVector dl(p1); dl -= p0;
    length = dl.GetMagnitude();
    ratio = std::sqrt(eigenValues.GetY() / eigenValues.GetX());
    transverse = eigenVectors[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeAssociationAlgorithm::MakeCone(const CartesianVector &start, const CartesianVector &finish, const CartesianVector &transverse,
    const float ratio, CartesianVector &origin, CartesianVector &vertex1, CartesianVector &vertex2) const
{
    // if m_runCount is 1, build an extended cone, projected 50% beyond the end of the axis, if m_runCount is 2
    origin = start;
    CartesianVector dl(finish); dl -= start;
    const float length{dl.GetMagnitude()};
    CartesianVector pcaExt(finish);

    if (m_runCount == 0)
    {   // Build an extended cone, projected 50% beyond the end of the axis
        dl *= 0.5f;
        pcaExt += dl;
        vertex1 = transverse;
        vertex1 *= 1.5f * length * ratio;
    }
    else
    {   // Build an extended cone, increasing the transverse extent by a factor of 2
        vertex1 = transverse;
        vertex1 *= length * ratio * 2;
    }
    vertex2 -= vertex1;
    vertex1 += pcaExt; vertex2 += pcaExt;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeAssociationAlgorithm::Run()
{
    m_runCount = 0;
    ClusterList allClusters;
    const ClusterList *pClusterList{nullptr};

    for (const std::string listName : m_inputListNames)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName,
            pClusterList));
        if (pClusterList && !pClusterList->empty())
        {
            m_viewToListMap[LArClusterHelper::GetClusterHitType(pClusterList->front())] = listName;
            for (const Cluster *const pCluster : *pClusterList)
                allClusters.emplace_back(pCluster);
        }
    }
    if (allClusters.empty())
    {
        std::cout << "ConeAssoicationAlgorithm: Error - no clusters found" << std::endl;
        return STATUS_CODE_NOT_INITIALIZED;
    }

    int repeat{0};
    while (true)
    {
        m_clusterToSortedHitsMap.clear();
        m_clustersU.clear(); m_clustersV.clear(); m_clustersW.clear();
        for (const Cluster *pCluster : allClusters)
        {
            const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
            if (view == TPC_VIEW_U)
                m_clustersU.emplace_back(pCluster);
            else if (view == TPC_VIEW_V)
                m_clustersV.emplace_back(pCluster);
            else if (view == TPC_VIEW_W)
                m_clustersW.emplace_back(pCluster);

            // TODO: May want to consider expanding this map to all clusters, not just clean clusters
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
            std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPositionInX);
            m_clusterToSortedHitsMap[pCluster] = caloHitVector;
        }

        // Collect the clean clusters in all views and sort them - doesn't matter that views are mixed here
        ClusterVector unsortedVector, clusterVector;
        this->GetListOfCleanClusters(&allClusters, unsortedVector);
        this->GetSortedListOfCleanClusters(unsortedVector, clusterVector);

        ClusterMergeMap clusterMergeMap;

        std::cout << "Pre #: " << m_runCount << " Repeat: " << repeat << std::endl;
        this->PopulateClusterMergeMap(clusterVector, clusterMergeMap);
        std::cout << "Post #: " << m_runCount << " Merged: " << !clusterMergeMap.empty() << " Repeat: " << repeat <<  std::endl;

        if (clusterMergeMap.empty() && m_runCount > 0)
            break;

        if (repeat < 3 && m_runCount > 1)
        {
            ++repeat;
            m_runCount = 0;
        }
        else if (repeat == 3 && m_runCount > 1)
        {
            break;
        }

        // Here we care about the views and Need to separate them
        std::map<HitType, ClusterVector> clusterVectors;
        std::map<HitType, ClusterMergeMap> clusterMergeMaps;
        for (const Cluster *pCluster : clusterVector)
            clusterVectors[LArClusterHelper::GetClusterHitType(pCluster)].emplace_back(pCluster);
        for (auto [ pSeed, targets ] : clusterMergeMap)
        {
            const HitType view{LArClusterHelper::GetClusterHitType(pSeed)};
            clusterMergeMaps[view][pSeed] = targets;
        }
        for (auto [ view, map ] : clusterMergeMaps)
        {
            try
            {
                // This is needed because the parent class makes assumptions about the list it is working with
                m_inputClusterListName = m_viewToListMap[view];
                if (!map.empty())
                    this->MergeClusters(clusterVectors[view], map);
            }
            catch (...)
            {
                std::cout << "Boom!" << std::endl;
            }
        }

        // Repopulate the cluster list
        allClusters.clear();
        for (const std::string listName : m_inputListNames)
        {
            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName,
                pClusterList));
            if (pClusterList && !pClusterList->empty())
            {
                m_viewToListMap[LArClusterHelper::GetClusterHitType(pClusterList->front())] = listName;
                for (const Cluster *const pCluster : *pClusterList)
                    allClusters.emplace_back(pCluster);
            }
        }
        if (allClusters.empty())
        {
            std::cout << "ConeAssoicationAlgorithm: Error - no clusters found" << std::endl;
            return STATUS_CODE_NOT_INITIALIZED;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputListNames", m_inputListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers",
        m_minClusterLayers));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize",
        m_visualize));

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
