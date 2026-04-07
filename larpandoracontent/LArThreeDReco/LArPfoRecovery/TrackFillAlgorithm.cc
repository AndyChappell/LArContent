/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackFillAlgorithm.cc
 *
 *  @brief  Implementation of the particle recovery algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPlaneContextObject.h"

#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackFillAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

TrackFillAlgorithm::TrackFillAlgorithm() :
    m_inputCaloHitListName("CaloHitList2D")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackFillAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));
    ViewToKDTreeMap kdTreeMap;
    this->FillKDTrees(*pCaloHitList, kdTreeMap);

    // Loop over all clusters to get a map from hits to their existing cluster
    HitToClusterToMap hitToClusterMap;
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, clusterListName, pClusterList))
            continue;
        this->FillHitToClusterMap(*pClusterList, hitToClusterMap);
    }

    // Get the list of track-like PFOs
    const PfoList *pPfoList(nullptr);
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList))
        return STATUS_CODE_SUCCESS;
    this->FillTracks(*pPfoList, hitToClusterMap, kdTreeMap);


    // Loop over the PFOs and project 3D hits from two views into the third view to look for unassociated hits
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    //    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsU, "U", RED));
    //        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "Merge U", RED, 2.f));
    //    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackFillAlgorithm::FillHitToClusterMap(const ClusterList &clusterList, HitToClusterToMap &hitToClusterMap) const
{
    for (const Cluster *const pCluster : clusterList)
    {
        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);
        for (const CaloHit *const pCaloHit : clusterHits)
            hitToClusterMap[pCaloHit] = pCluster;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackFillAlgorithm::FitAndOrderCluster(const Cluster *const pCluster, CaloHitList &orderedHits) const
{
    if (pCluster->GetNCaloHits() < 3)
    {
        LArClusterHelper::GetAllHits(pCluster, orderedHits);
        return false;
    }
    const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
    TwoDSlidingFitResult sfr(pCluster, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
    LArClusterHelper::OrderHitsAlongTrajectory(pCluster, sfr, orderedHits);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackFillAlgorithm::FillTracks(const PfoList &pfoList, HitToClusterToMap &hitToClusterMap, ViewToKDTreeMap &kdTreeMap) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    std::unordered_map<const CaloHit *, int> prevMergesMap;
    bool mergesMade{false};
    do
    {
        mergesMade = false;
        std::unordered_map<const CaloHit *, int> mergesMap;
        for (const Pfo *const pPfo : pfoList)
        {
            for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                ClusterList pfoClusters;
                LArPfoHelper::GetClusters(pPfo, view, pfoClusters);
                if (!pfoClusters.empty())
                {
                    const Cluster *pCluster{pfoClusters.front()};
                    CaloHitList clusterHits;
                    if (this->FitAndOrderCluster(pCluster, clusterHits))
                    {
                        CaloHitList mergeHits;
                        for (auto iter = clusterHits.begin(); iter != std::prev(clusterHits.end()); ++iter)
                        {
                            const CaloHit *const pCaloHit1{*iter}, *const pCaloHit2{*(std::next(iter))};
                            const CartesianVector position1{pCaloHit1->GetPositionVector()}, position2{pCaloHit2->GetPositionVector()};
                            const CartesianVector midPoint{(position1 + position2) * 0.5f};
                            const float dx{std::abs(position1.GetX() - position2.GetX())}, dz{std::abs(position1.GetZ() - position2.GetZ())};
                            KDTreeBox searchHits(build_2d_kd_search_region(midPoint, dx, dz));
                            NodeList found;
                            kdTreeMap.at(view).search(searchHits, found);
                            if (found.size() > 2)
                            {
                                for (const auto &hit : found)
                                {
                                    const CaloHit *const pCaloHit(hit.data);
                                    if (pCaloHit == pCaloHit1 || pCaloHit == pCaloHit2)
                                        continue;
                                    if (LArClusterHelper::HasBlockedPath({pCaloHit}, pCaloHit1, pCaloHit2))
                                        mergeHits.emplace_back(pCaloHit);
                                }
                            }
                        }
                        if (!mergeHits.empty())
                        {
                            for (const CaloHit *const pCaloHit : mergeHits)
                            {
                                const Cluster *pOldCluster{hitToClusterMap.at(pCaloHit)};
                                PandoraContentApi::RemoveFromCluster(*this, pOldCluster, pCaloHit);
                                PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit);
                                hitToClusterMap[pCaloHit] = pCluster;

                                if (mergesMap.find(pCaloHit) == mergesMap.end())
                                    mergesMap[pCaloHit] = 0;
                                mergesMap[pCaloHit]++;
                            }
                            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &clusterHits, "original", BLUE));
                            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mergeHits, "new", RED));
                            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                        }
                    }
                }
            }
        }
        // Look for unambiguous merges (i.e. hits that only get moved once). We want to repeat this process until no unambiguous merges remain
        for (const auto &[pCaloHit, numMerges] : mergesMap)
        {
            if (numMerges && prevMergesMap.find(pCaloHit) == prevMergesMap.end())
            {
                mergesMade = true;
                std::cout << "Found a new merge" << std::endl;
                break;
            }
        }
        prevMergesMap = std::move(mergesMap);
    }
    while (mergesMade);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackFillAlgorithm::FillKDTrees(const CaloHitList &caloHitList, ViewToKDTreeMap &kdTreeMap) const
{
    std::unordered_map<const HitType, CaloHitList> caloHitMap;
    for (const CaloHit *const pCaloHit : caloHitList)
        caloHitMap[pCaloHit->GetHitType()].emplace_back(pCaloHit);

    for (HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        NodeList nodeList;
        KDTreeBox bounds{fill_and_bound_2d_kd_tree(caloHitMap[view], nodeList)};
        kdTreeMap[view].build(nodeList, bounds);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackFillAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName",
        m_inputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
