/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/EnvelopeAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the proximity association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/EnvelopeAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EnvelopeAssociationAlgorithm::AssociationCandidate::AssociationCandidate(const Cone &cone, const ClusterList &clusters) :
    m_cone(cone),
    m_clusters(clusters)
{
    for (const Cluster *pCluster : clusters)
    {
        const OrderedCaloHitList &orderedCaloHits(pCluster->GetOrderedCaloHitList());
        orderedCaloHits.FillCaloHitList(m_hits);
    }
    m_sortedHits = CaloHitVector(m_hits.begin(), m_hits.end());
    std::sort(m_sortedHits.begin(), m_sortedHits.end(), LArClusterHelper::SortHitsByPositionInX);
    m_xMin = m_sortedHits.front()->GetPositionVector().GetX();
    m_xMax = m_sortedHits.back()->GetPositionVector().GetX();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnvelopeAssociationAlgorithm::AssociationCandidate::GetOverlapFraction(const AssociationCandidate &other, float &minOverlapFraction,
    float &maxOverlapFraction) const
{
    const float thisRange{m_xMax - m_xMin};
    const float otherRange{other.m_xMax - other.m_xMin};
    const float overlapStart{std::max(m_xMin, other.m_xMin)};
    const float overlapFinish{std::min(m_xMax, other.m_xMax)};
    if (overlapFinish >= overlapStart)
    {
        const float overlapSize{overlapFinish - overlapStart};
        if (overlapSize > std::numeric_limits<float>::epsilon())
        {
            // ATTN - Neither range can be zero if we're here, so need to check
            const float overlap1{overlapSize / thisRange}, overlap2{overlapSize / otherRange};
            minOverlapFraction = std::min(overlap1, overlap2);
            maxOverlapFraction = std::max(overlap1, overlap2);
        }
        else
        {   // At least one candidate is isochronous
            if (thisRange < std::numeric_limits<float>::epsilon() && otherRange < std::numeric_limits<float>::epsilon())
            {   // Both isochronous
                minOverlapFraction = 1.f;
                maxOverlapFraction = 1.f;
            }
            else
            {   // One candidate is isochronous
                minOverlapFraction = 0.f;
                maxOverlapFraction = 1.f;
            }
        }
    }
    else
    {   // No overlap at all
        minOverlapFraction = 0.f;
        maxOverlapFraction = 0.f;
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnvelopeAssociationAlgorithm::OverlapCandidates::OverlapCandidates(const Pandora &pandora, const AssociationCandidate &candidate1,
    const AssociationCandidate &candidate2) :
    m_nViews(2),
    m_candidateClusters1(candidate1.GetClusters()),
    m_candidateClusters2(candidate2.GetClusters()),
    m_chi2(std::numeric_limits<float>::max())
{
    const CaloHitVector &hits1{candidate1.GetSortedCaloHits()};
    const CaloHitVector &hits2{candidate2.GetSortedCaloHits()};

    const float overlapStart{std::max(candidate1.GetMinX(), candidate2.GetMinX())};
    const float overlapFinish{std::min(candidate1.GetMaxX(), candidate2.GetMaxX())};
    // Bins are 0.1 cm wide - ensure 1 bin in case start == finish (isochronous case)
    const int N{std::max(1, static_cast<int>(std::ceil(10.f * (overlapFinish - overlapStart))))};
    CaloHitList *pBinnedHits1{new CaloHitList[N]}, *pBinnedHits2{new CaloHitList[N]};

    for (const CaloHit *pCaloHit : hits1)
    {
        const float hitX{pCaloHit->GetPositionVector().GetX()};
        const int bin{static_cast<int>(std::floor(10.f * (hitX - overlapStart)))};
        if (bin >= 0 && bin < N)
            pBinnedHits1[bin].emplace_back(pCaloHit);
        // ATTN - assumes hits are sorted
        if (bin >= N)
            break;
    }

    for (const CaloHit *pCaloHit : hits2)
    {
        const float hitX{pCaloHit->GetPositionVector().GetX()};
        const int bin{static_cast<int>(std::floor(10.f * (hitX - overlapStart)))};
        if (bin >= 0 && bin < N)
            pBinnedHits2[bin].emplace_back(pCaloHit);
        // ATTN - assumes hits are sorted
        if (bin >= N)
            break;
    }

    int goodBins{0};
    int usedBins{0};
    CartesianVector dummy(0.f, 0.f, 0.f);
    for (int i = 0; i < N; ++i)
    {
        if (!pBinnedHits1[i].empty() && !pBinnedHits2[i].empty())
        {
            float minChi2{std::numeric_limits<float>::max()};
            for (const CaloHit *pCaloHit1 : pBinnedHits1[i])
            {
                for (const CaloHit *pCaloHit2 : pBinnedHits2[i])
                {
                    float chi2{0.f};
                    LArGeometryHelper::MergeTwoPositions3D(pandora, pCaloHit1->GetHitType(), pCaloHit2->GetHitType(),
                        pCaloHit1->GetPositionVector(), pCaloHit2->GetPositionVector(), dummy, chi2);
                    if (chi2 < minChi2)
                        minChi2 = chi2;
                }
            }
            if (minChi2 < std::numeric_limits<float>::max())
            {
                std::cout << i << " of " << N << ": " << minChi2 << std::endl;
                if (minChi2 < 0.1f)
                    ++goodBins;
                ++usedBins;
            }
        }
    }
    m_chi2 = usedBins > 0 ? goodBins / static_cast<float>(usedBins) : 0.f;

    // Note, may not want to delete at this stage in the future, might defer to destructor
    delete[] pBinnedHits1;
    delete[] pBinnedHits2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

EnvelopeAssociationAlgorithm::OverlapCandidates::OverlapCandidates(const Pandora &pandora, const AssociationCandidate &candidate1,
    const AssociationCandidate &candidate2, const AssociationCandidate &candidate3) :
    m_nViews(3),
    m_candidateClusters1(candidate1.GetClusters()),
    m_candidateClusters2(candidate2.GetClusters()),
    m_candidateClusters3(candidate3.GetClusters()),
    m_chi2(std::numeric_limits<float>::max())
{
    const CaloHitVector &hits1{candidate1.GetSortedCaloHits()};
    const CaloHitVector &hits2{candidate2.GetSortedCaloHits()};
    const CaloHitVector &hits3{candidate3.GetSortedCaloHits()};

    const float overlapStart{std::max({candidate1.GetMinX(), candidate2.GetMinX(), candidate3.GetMinX()})};
    const float overlapFinish{std::min({candidate1.GetMaxX(), candidate2.GetMaxX(), candidate3.GetMaxX()})};
    // Bins are 0.1 cm wide - ensure 1 bin in case start == finish (isochronous case)
    const int N{std::max(1, static_cast<int>(std::ceil(10.f * (overlapFinish - overlapStart))))};
    CaloHitList *pBinnedHits1{new CaloHitList[N]}, *pBinnedHits2{new CaloHitList[N]}, *pBinnedHits3{new CaloHitList[N]};

    for (const CaloHit *pCaloHit : hits1)
    {
        const float hitX{pCaloHit->GetPositionVector().GetX()};
        const int bin{static_cast<int>(std::floor(10.f * (hitX - overlapStart)))};
        if (bin >= 0 && bin < N)
            pBinnedHits1[bin].emplace_back(pCaloHit);
        // ATTN - assumes hits are sorted
        if (bin >= N)
            break;
    }

    for (const CaloHit *pCaloHit : hits2)
    {
        const float hitX{pCaloHit->GetPositionVector().GetX()};
        const int bin{static_cast<int>(std::floor(10.f * (hitX - overlapStart)))};
        if (bin >= 0 && bin < N)
            pBinnedHits2[bin].emplace_back(pCaloHit);
        // ATTN - assumes hits are sorted
        if (bin >= N)
            break;
    }

    for (const CaloHit *pCaloHit : hits3)
    {
        const float hitX{pCaloHit->GetPositionVector().GetX()};
        const int bin{static_cast<int>(std::floor(10.f * (hitX - overlapStart)))};
        if (bin >= 0 && bin < N)
            pBinnedHits3[bin].emplace_back(pCaloHit);
        // ATTN - assumes hits are sorted
        if (bin >= N)
            break;
    }

    int goodBins{0};
    int usedBins{0};
    CartesianVector dummy(0.f, 0.f, 0.f);
    for (int i = 0; i < N; ++i)
    {
        if (!pBinnedHits1[i].empty() && !pBinnedHits2[i].empty() && !pBinnedHits3[i].empty())
        {
            float minChi2{std::numeric_limits<float>::max()};
            for (const CaloHit *pCaloHit1 : pBinnedHits1[i])
            {
                for (const CaloHit *pCaloHit2 : pBinnedHits2[i])
                {
                    for (const CaloHit *pCaloHit3 : pBinnedHits3[i])
                    {
                        float chi2{0.f};
                        LArGeometryHelper::MergeThreePositions3D(pandora, pCaloHit1->GetHitType(), pCaloHit2->GetHitType(),
                            pCaloHit3->GetHitType(), pCaloHit1->GetPositionVector(), pCaloHit2->GetPositionVector(),
                            pCaloHit3->GetPositionVector(), dummy, chi2);
                        if (chi2 < minChi2)
                            minChi2 = chi2;
                    }
                }
            }
            if (minChi2 < std::numeric_limits<float>::max())
            {
                if (minChi2 < 0.1f)
                    ++goodBins;
                ++usedBins;
            }
        }
        // Might want to consider and else block with 2 view 3D position determination where 3 views aren't available
    }
    m_chi2 = usedBins > 0 ? goodBins / static_cast<float>(usedBins) : 0.f;

    // Note, may not want to delete at this stage in the future, might defer to destructor
    delete[] pBinnedHits1;
    delete[] pBinnedHits2;
    delete[] pBinnedHits3;
}

//------------------------------------------------------------------------------------------------------------------------------------------


void EnvelopeAssociationAlgorithm::OverlapCandidates::AddToMergeMap(ClusterMergeMap &clusterMergeMap) const
{
    const Cluster *pSeed1{m_candidateClusters1.front()};
    for (const Cluster *pCluster : m_candidateClusters1)
    {
        if (pCluster != pSeed1 && clusterMergeMap.find(pCluster) == clusterMergeMap.end())
        {
            clusterMergeMap[pSeed1].emplace_back(pCluster);
            clusterMergeMap[pCluster].emplace_back(pSeed1);
        }
    }
    const Cluster *pSeed2{m_candidateClusters2.front()};
    for (const Cluster *pCluster : m_candidateClusters2)
    {
        if (pCluster != pSeed2 && clusterMergeMap.find(pCluster) == clusterMergeMap.end())
        {
            clusterMergeMap[pSeed2].emplace_back(pCluster);
            clusterMergeMap[pCluster].emplace_back(pSeed2);
        }
    }

    if (m_nViews == 3)
    {
        const Cluster *pSeed3{m_candidateClusters3.front()};
        for (const Cluster *pCluster : m_candidateClusters3)
        {
            if (pCluster != pSeed3 && clusterMergeMap.find(pCluster) == clusterMergeMap.end())
            {
                clusterMergeMap[pSeed3].emplace_back(pCluster);
                clusterMergeMap[pCluster].emplace_back(pSeed3);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnvelopeAssociationAlgorithm::OverlapCandidates::SharesSeed(const OverlapCandidates &other) const
{
    const Cluster *pThisSeed1{m_candidateClusters1.front()}, *pThisSeed2{m_candidateClusters2.front()};
    const Cluster *pOtherSeed1{other.m_candidateClusters1.front()}, *pOtherSeed2{other.m_candidateClusters2.front()};
    if (m_nViews == 2)
    {
        return pThisSeed1 == pOtherSeed1 || pThisSeed1 == pOtherSeed2 || pThisSeed2 == pOtherSeed1 || pThisSeed2 == pOtherSeed2;
    }
    else
    {
        const Cluster *pThisSeed3{m_candidateClusters3.front()};
        const Cluster *pOtherSeed3{other.m_candidateClusters3.front()};
        return pThisSeed1 == pOtherSeed1 || pThisSeed1 == pOtherSeed2 || pThisSeed1 == pOtherSeed3 ||
            pThisSeed2 == pOtherSeed1 || pThisSeed2 == pOtherSeed2 || pThisSeed2 == pOtherSeed3 ||
            pThisSeed3 == pOtherSeed1 || pThisSeed3 == pOtherSeed2 || pThisSeed3 == pOtherSeed3;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnvelopeAssociationAlgorithm::EnvelopeAssociationAlgorithm() :
    m_minClusterLayers{1},
    m_minSeedCaloHits{30},
    m_runCount{0},
    m_visualize{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnvelopeAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
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

void EnvelopeAssociationAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    if (m_runCount > 1)
        return;
    if (m_runCount == 0)
        this->AssociateClusters(clusterVector, clusterMergeMap, m_minSeedCaloHits);
    else
        this->AssociateClusters(clusterVector, clusterMergeMap, 3 * m_minSeedCaloHits);
    ++m_runCount;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnvelopeAssociationAlgorithm::AssociateClusters(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap,
    const unsigned int minSeedCaloHits) const
{
    // Separate the clusters into their respective views
    std::map<HitType, ClusterList> viewToClustersMap;
    for (const Cluster *pCluster : clusterVector)
        viewToClustersMap[LArClusterHelper::GetClusterHitType(pCluster)].emplace_back(pCluster);
    // Might want to throw here
    if (viewToClustersMap.empty())
        return;

    std::list<HitType> availableViews;
    for (const auto [ view, list ] : viewToClustersMap)
    {   (void)list;
        availableViews.emplace_back(view);
    }
    ViewToAssociationCandidatesMap viewToCandidatesMap;
    for (const auto [ currentView, currentList ] : viewToClustersMap)
    {
        ClusterList seedClusters, targetClusters;
        for (const Cluster *pCluster : currentList)
        {
            if (pCluster->GetNCaloHits() > minSeedCaloHits)
                seedClusters.emplace_back(pCluster);
            else
                targetClusters.emplace_back(pCluster);
        }
        for (const Cluster *pSeed : seedClusters)
        {
            Cone cone(this->GetBoundingCone(pSeed));
            ClusterList containedClusters{pSeed};
            for (const Cluster * pTarget : targetClusters)
            {
                if (!this->IsClusterContained(cone, pTarget))
                    continue;
                containedClusters.emplace_back(pTarget);
            }
            AssociationCandidate candidate(cone, containedClusters);
            viewToCandidatesMap[currentView].emplace_back(candidate);
        }
    }

    // Collect up the overlap candidates - can have competing views at this stage
    std::list<OverlapCandidates> overlapCandidatesList;
    float minOverlapVW{0.f}, maxOverlapVW{0.f};
    float minOverlapUW{0.f}, maxOverlapUW{0.f};
    float minOverlapUV{0.f}, maxOverlapUV{0.f};
    for (const AssociationCandidate &candidateW : viewToCandidatesMap[TPC_VIEW_W])
    {
        for (const AssociationCandidate &candidateV : viewToCandidatesMap[TPC_VIEW_V])
        {
            candidateW.GetOverlapFraction(candidateV, minOverlapVW, maxOverlapVW);
            for (const AssociationCandidate &candidateU : viewToCandidatesMap[TPC_VIEW_U])
            {
                candidateV.GetOverlapFraction(candidateU, minOverlapUV, maxOverlapUV);
                candidateW.GetOverlapFraction(candidateU, minOverlapUW, maxOverlapUW);
                const float minOverlap{std::min({minOverlapVW, minOverlapUV, minOverlapUW})};
                const float maxOverlap{std::max({maxOverlapVW, maxOverlapUV, maxOverlapUW})};
                bool overlaps{maxOverlap > 0.66f && minOverlap > 0.5f};
                if (!overlaps && maxOverlap == 1.f)
                {   // There's at least one isochronous case, all 3 candidates need maxOverlap == 1
                    overlaps = maxOverlapUV == 1.f && maxOverlapUW == 1.f && maxOverlapVW == 1.f;
                }
                if (overlaps)
                {
                    OverlapCandidates overlap(this->GetPandora(), candidateW, candidateV, candidateU);
                    if (overlap.GetChiSquared() > 0.66f)
                    {
                        std::cout << "Matched " << overlap.GetChiSquared() << std::endl;
                        overlapCandidatesList.emplace_back(overlap);
/*                        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidateU.GetClusters(), "1", RED, false));
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidateV.GetClusters(), "2", GREEN, false));
                        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidateW.GetClusters(), "3", BLUE, false));
                        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
                    }
                }
            }
        }
    }

    // Try two view matches
    for (const AssociationCandidate &candidateW : viewToCandidatesMap[TPC_VIEW_W])
    {
        for (const AssociationCandidate &candidateV : viewToCandidatesMap[TPC_VIEW_V])
        {
            float minOverlap{0.f}, maxOverlap{0.f};
            candidateW.GetOverlapFraction(candidateV, minOverlap, maxOverlap);
            const bool overlaps{(maxOverlap > 0.66f && minOverlap > 0.5f) || maxOverlap == 1.f};
            if (overlaps)
            {
                OverlapCandidates overlap(this->GetPandora(), candidateW, candidateV);
                if (overlap.GetChiSquared() > 0.66f)
                {
                    std::cout << "Matched two view " << overlap.GetChiSquared() << std::endl;
                    overlapCandidatesList.emplace_back(overlap);
/*                    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
                    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidateV.GetClusters(), "V", GREEN, false));
                    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidateW.GetClusters(), "W", BLUE, false));
                    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
                }
            }
        }
    }

    overlapCandidatesList.sort();
    std::list<OverlapCandidates> mergeList;
    for (const OverlapCandidates &candidates : overlapCandidatesList)
    {
        bool isGood{true};
        for (const OverlapCandidates &goodCandidates : mergeList)
        {
            if (candidates.SharesSeed(goodCandidates))
            {
                isGood = false;
                break;
            }
        }
        if (isGood)
        {
            mergeList.emplace_back(candidates);
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidates.m_candidateClusters1, "1", RED, false));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidates.m_candidateClusters2, "2", GREEN, false));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &candidates.m_candidateClusters3, "3", BLUE, false));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
    // Flag the best candidates for merging
    for (const OverlapCandidates &candidates : mergeList)
        candidates.AddToMergeMap(clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnvelopeAssociationAlgorithm::IsClusterContained(const Cone &cone, const Cluster *const pCluster) const
{
    CartesianVector ab(cone.m_b); ab -= cone.m_a;
    CartesianVector ac(cone.m_c); ac -= cone.m_a;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);
    CaloHitList containedHits, uncontainedHits;
    for (const CaloHit *pCaloHit : caloHits)
    {
        CartesianVector ap(pCaloHit->GetPositionVector()); ap -= cone.m_a;
        if (this->IsPointContained(ab, ac, ap))
            containedHits.emplace_back(pCaloHit);
        else
            uncontainedHits.emplace_back(pCaloHit);
    }

    return !containedHits.empty() && containedHits.size() >= uncontainedHits.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnvelopeAssociationAlgorithm::IsPointContained(const CartesianVector &ab, const CartesianVector &ac, const CartesianVector &ap) const
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

EnvelopeAssociationAlgorithm::Cone EnvelopeAssociationAlgorithm::GetBoundingCone(const Cluster *const pCluster) const
{
    // Retrieve the hits for PCA
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);

    // Get the eigen vectors for this cluster
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(caloHits, centroid, eigenValues, eigenVectors);

    // Project the extremal hits of the cluster onto the principal axis
    CartesianVector p0(0.f, 0.f, 0.f), p1(0.f, 0.f, 0.f);
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
    const float length{dl.GetMagnitude()};

    const float eigenRatio{std::sqrt(eigenValues.GetY() / eigenValues.GetX())};
    // Define a cone that (approximately) envelopes the cluster
    CartesianVector p1r(eigenVectors[1]), p1l(0.f, 0.f, 0.f);
    p1r *= length * eigenRatio;
    p1l -= p1r;
    p1r += p1; p1l += p1;

    CartesianVector p2(p1);
    CartesianVector p2r(eigenVectors[1]), p2l(0.f, 0.f, 0.f);
    if (m_runCount > 0)
    {   // Broader cone surrounding the seed cluster
        p2r *= 3.5f * length * eigenRatio;
        p2l -= p2r;
        p2r += p2; p2l += p2;
    }
    else
    {
        // Define an extended cone beyond the end of the seed cluster
        p2 += dl;
        p2r *= 2 * length * eigenRatio;
        p2l -= p2r;
        p2r += p2; p2l += p2;
/*        p2r *= 3.5f * length * eigenRatio;
        p2l -= p2r;
        p2r += p2; p2l += p2;*/
    }

    if (m_visualize && LArClusterHelper::GetClusterHitType(pCluster) == HitType::TPC_VIEW_W)
    {
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p1r, "A", BLUE, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p1l, "B", BLUE, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p0, "C", BLUE, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p2r, "D", RED, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p2l, "E", RED, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p2r, &p2l, "H", RED, 3, 1));
        PANDORA_MONITORING_API(Pause(this->GetPandora()));
    }

    if (m_visualize && LArClusterHelper::GetClusterHitType(pCluster) == HitType::TPC_VIEW_W)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "Seed", BLUE));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p2l, "AB", BLACK, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p2r, "AC", BLACK, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p2l, &p2r, "BC", BLACK, 3, 1));
        PANDORA_MONITORING_API(Pause(this->GetPandora()));
    }

    Cone cone(p0, p2l, p2r);    // origin, 'left' vertex, 'right' vertex
    return cone;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnvelopeAssociationAlgorithm::Run()
{
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
        std::cout << "EnvelopeAssoicationAlgorithm: Error - no clusters found" << std::endl;
        return STATUS_CODE_NOT_INITIALIZED;
    }

    while (true)
    {
        // Collect the clean clusters in all views and sort them - doesn't matter that views are mixed here
        ClusterVector unsortedVector, clusterVector;
        this->GetListOfCleanClusters(&allClusters, unsortedVector);
        this->GetSortedListOfCleanClusters(unsortedVector, clusterVector);

        ClusterMergeMap clusterMergeMap;
        this->PopulateClusterMergeMap(clusterVector, clusterMergeMap);

        if (clusterMergeMap.empty())
            break;

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
            std::cout << "EnvelopeAssoicationAlgorithm: Error - no clusters found" << std::endl;
            return STATUS_CODE_NOT_INITIALIZED;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnvelopeAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputListNames", m_inputListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers",
        m_minClusterLayers));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSeedCaloHits",
        m_minSeedCaloHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize",
        m_visualize));

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
