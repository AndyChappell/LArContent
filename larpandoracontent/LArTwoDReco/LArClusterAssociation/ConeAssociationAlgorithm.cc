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

void ConeAssociationAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    if (m_runCount > 1)
        return;
    // Separate the clusters into their respective views
    ClusterList clusterListU, clusterListV, clusterListW;
    for (const Cluster *pCluster : clusterVector)
    {
        HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        if (view == TPC_VIEW_U)
            clusterListU.emplace_back(pCluster);
        else if (view == TPC_VIEW_V)
            clusterListV.emplace_back(pCluster);
        else if (view == TPC_VIEW_W)
            clusterListW.emplace_back(pCluster);
    }

    if (m_runCount == 0)
    {   // Iteration 1 uses x-overlap information to consolidate clusters ahead of PCA step
        if (!(clusterListU.empty() || clusterListV.empty() || clusterListW.empty()))
        {
            ViewClusterVector viewClusterVector;
            this->AssociateClusters(clusterListU, clusterListV, clusterListW, viewClusterVector);
            this->ConsolidateClusters(viewClusterVector, clusterMergeMap);
            for (const ViewCluster *pViewCluster : viewClusterVector)
                delete pViewCluster;
        }
    }
    else if (m_runCount == 1)
    {   // Iteration 2 builds PCA cones to grow clusters
        if (!(clusterListU.empty() || clusterListV.empty() || clusterListW.empty()))
        {
            ViewClusterVector viewClusterVector;
            this->AssociateClusters(clusterListU, clusterListV, clusterListW, viewClusterVector);
            this->GrowClusters(viewClusterVector, clusterMergeMap);
            for (const ViewCluster *pViewCluster : viewClusterVector)
                delete pViewCluster;
        }
    }
    ++m_runCount;
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

        std::cout << "Start 1: (" << pca1Start.GetX() << "," << pca1Start.GetZ() << ") 2: (" << pca2Start.GetX() << "," << pca2Start.GetZ() <<
            ")" << std::endl;
        std::cout << "Finish 1: (" << pca1Finish.GetX() << "," << pca1Finish.GetZ() << ") 2: (" << pca2Finish.GetX() << "," << pca2Finish.GetZ() <<
            ")" << std::endl;

        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca1Start, &pca1Finish, "1", ORANGE, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca2Start, &pca2Finish, "2", MAGENTA, 3, 1));

        if (pca1Start.GetX() != pca1Finish.GetX() && pca2Start.GetX() != pca2Finish.GetX())
        {   // Neither axis is isochronous, endpoints can be matched on x alone
            if (pca1Start.GetX() == pca2Start.GetX())
            {   // Matching start to start and finish to finish
                if (views[0] == TPC_VIEW_U)
                {
                    if (views[1] == TPC_VIEW_V)
                    {   // UV inputs
                        const float us{pca1Start.GetZ()}, vs{pca2Start.GetZ()};
                        const float ws{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(us, vs))};
                        CartesianVector pca3Start(pca1Start.GetX(), 0.f, ws);
                        const float uf{pca1Finish.GetZ()}, vf{pca2Finish.GetZ()};
                        const float wf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(uf, vf))};
                        CartesianVector pca3Finish(pca1Finish.GetX(), 0.f, wf);
                        //CartesianVector pca3DStart(0.f, 0.f, 0.f);
                        //float chi2{std::numeric_limits::max()}
                        //LArGeometryHelper::MergeThreePositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, pca1Start, pca2Start, pca3Start,
                        //    pca3DStart, chi2);
                        //std::cout << "Proj W chi2: " << chi2 << std::endl;
                        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjW", BLACK, 3, 1));

                        CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                        CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                        const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                        std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                        if (cosOpenAngle > 0.98f)
                        {
                            std::cout << "Worth pursuing" << std::endl;
                        }
                    }
                    else
                    {   // UW inputs
                        const float us{pca1Start.GetZ()}, ws{pca2Start.GetZ()};
                        const float vs{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(ws, us))};
                        CartesianVector pca3Start(pca1Start.GetX(), 0.f, vs);
                        const float uf{pca1Finish.GetZ()}, wf{pca2Finish.GetZ()};
                        const float vf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(wf, uf))};
                        CartesianVector pca3Finish(pca1Finish.GetX(), 0.f, vf);
                        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjV", BLACK, 3, 1));

                        CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                        CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                        const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                        std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                        if (cosOpenAngle > 0.98f)
                        {
                            std::cout << "Worth pursuing" << std::endl;
                        }
                    }
                }
                else
                {   // VW inputs
                    const float vs{pca1Start.GetZ()}, ws{pca2Start.GetZ()};
                    const float us{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vs, ws))};
                    CartesianVector pca3Start(pca1Start.GetX(), 0.f, us);
                    const float vf{pca1Finish.GetZ()}, wf{pca2Finish.GetZ()};
                    const float uf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vf, wf))};
                    CartesianVector pca3Finish(pca1Finish.GetX(), 0.f, uf);
                    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjU", BLACK, 3, 1));

                    CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                    CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                    const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                    std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                    if (cosOpenAngle > 0.98f)
                    {
                        std::cout << "Worth pursuing" << std::endl;
                    }
                }
            }
            else
            {   // Matching start1 to finish2 and start2 to finish1
                if (views[0] == TPC_VIEW_U)
                {
                    if (views[1] == TPC_VIEW_V)
                    {   // UV inputs
                        const float us{pca1Start.GetZ()}, vs{pca2Finish.GetZ()};
                        const float ws{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(us, vs))};
                        CartesianVector pca3Start(pca1Start.GetX(), 0.f, ws);
                        const float uf{pca1Finish.GetZ()}, vf{pca2Start.GetZ()};
                        const float wf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(uf, vf))};
                        CartesianVector pca3Finish(pca1Finish.GetX(), 0.f, wf);
                        //CartesianVector pca3DStart(0.f, 0.f, 0.f);
                        //float chi2{std::numeric_limits::max()}
                        //LArGeometryHelper::MergeThreePositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, pca1Start, pca2Start, pca3Start,
                        //    pca3DStart, chi2);
                        //std::cout << "Proj W chi2: " << chi2 << std::endl;
                        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjW", BLACK, 3, 1));

                        CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                        CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                        const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                        std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                        if (cosOpenAngle > 0.98f)
                        {
                            std::cout << "Worth pursuing" << std::endl;
                        }
                    }
                    else
                    {   // UW inputs
                        const float us{pca1Start.GetZ()}, ws{pca2Finish.GetZ()};
                        const float vs{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(ws, us))};
                        CartesianVector pca3Start(pca1Start.GetX(), 0.f, vs);
                        const float uf{pca1Finish.GetZ()}, wf{pca2Start.GetZ()};
                        const float vf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(wf, uf))};
                        CartesianVector pca3Finish(pca1Finish.GetX(), 0.f, vf);
                        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjV", BLACK, 3, 1));

                        CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                        CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                        const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                        std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                        if (cosOpenAngle > 0.98f)
                        {
                            std::cout << "Worth pursuing" << std::endl;
                        }
                    }
                }
                else
                {   // VW inputs
                    const float vs{pca1Start.GetZ()}, ws{pca2Finish.GetZ()};
                    const float us{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vs, ws))};
                    CartesianVector pca3Start(pca1Start.GetX(), 0.f, us);
                    const float vf{pca1Finish.GetZ()}, wf{pca2Start.GetZ()};
                    const float uf{static_cast<float>(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(vf, wf))};
                    CartesianVector pca3Finish(pca1Finish.GetX(), 0.f, uf);
                    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pca3Start, &pca3Finish, "ProjU", BLACK, 3, 1));

                    CartesianVector axis1(pca3Finish); axis1 -= pca3Start;
                    CartesianVector axis2(pcaXFinish); axis2 -= pcaXStart;
                    const float cosOpenAngle{std::abs(axis1.GetCosOpeningAngle(axis2))};
                    std::cout << "cosOpenAngle: " << cosOpenAngle << std::endl;
                    if (cosOpenAngle > 0.98f)
                    {
                        std::cout << "Worth pursuing" << std::endl;
                    }
                }
            }
        }
        else
        {   // At least one axis is isochronous, need to check all possible endpoint matches
            std::cout << "Isochronous case" << std::endl;
        }


        // Axes have been correctly extended at this point, next need to take the end points (first match based on x, if isochronous pick the best
        // projection chi2 to move forward) and project into the remaining view. Might want to apply transverse coords here and then generate 3D hits
        // from the vertices of the cone, before projecting the 3D points back into the 3 views and finding clusters that are mostly contained and 
        // also produce good 3D points if they were to be merged

        // Define a cone that (approximately) envelopes the cluster
/*        CartesianVector pca1Start(0.f, 0.f, 0.f), pca1Finish(0.f, 0.f, 0.f), p1r(0.f, 0.f, 0.f), p1l(0.f, 0.f, 0.f), transverse(0.f, 0.f, 0.f);

        p1r *= length * ratio;
        p1l -= p1r;
        p1r += p1; p1l += p1;

        for (HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            const Cluster *pCluster{viewCluster.GetCluster(view)};
            // Project the extremal hits of the cluster onto the principal axis
            CartesianVector p0(0.f, 0.f, 0.f), p1(0.f, 0.f, 0.f), p1r(0.f, 0.f, 0.f), p1l(0.f, 0.f, 0.f), transverse(0.f, 0.f, 0.f);
            float length{0.f}, ratio{1.f};
            this->GetConeParameters(pCluster, p0, p1, transverse, length, ratio);
            // Define a cone that (approximately) envelopes the cluster
            p1r *= length * ratio;
            p1l -= p1r;
            p1r += p1; p1l += p1;

            //if (m_visualize && LArClusterHelper::GetClusterHitType(pCluster) == HitType::TPC_VIEW_W)
            {
                Color color{view == TPC_VIEW_W ? BLUE : view == TPC_VIEW_V ? GREEN : RED};
                CaloHitList &caloHits{view == TPC_VIEW_W ? caloHitsW : view == TPC_VIEW_V ? caloHitsV : caloHitsU};
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, std::to_string(view), GRAY));
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p1r, "A", color, 3, 1));
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p1l, "B", color, 3, 1));
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p0, "C", color, 3, 1));
            }
        }*/
        PANDORA_MONITORING_API(Pause(this->GetPandora()));


        /*
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
        }*/






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

StatusCode ConeAssociationAlgorithm::Run()
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
        std::cout << "ConeAssoicationAlgorithm: Error - no clusters found" << std::endl;
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
