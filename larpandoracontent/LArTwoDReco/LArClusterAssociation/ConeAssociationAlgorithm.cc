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
    m_minHitsInOverlapRegion{std::numeric_limits<int>::max()}
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
    if (m_runCount == 0)
    {
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
        if (!(clusterListU.empty() || clusterListV.empty() || clusterListW.empty()))
        {
            ViewClusterVector viewClusterVector;
            this->AssociateClusters(clusterListU, clusterListV, clusterListW, viewClusterVector);
            this->ConsolidateClusters(viewClusterVector, clusterMergeMap);
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
