/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.cc
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

#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

TrackRecoveryAlgorithm::TrackRecoveryAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::Run()
{
    // Get the list of track-like PFOs
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    // Get maps between PFOs, views, clusters and hits
    PfoToViewClusterMap pfoToViewClusterMap;
    ClusterToPfoMap clusterToPfoMap;
    ClusterToHitMap clusterToHitMap;
    ClusterToFitMap clusterToFitMap;
    HitToClusterToMap hitToClusterMap;
    ViewToHitsMap viewToHitsMap;
    for (const Pfo *const pPfo : *pPfoList)
    {
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, TPC_3D})
        {
            ClusterList pfoClusters;
            LArPfoHelper::GetClusters(pPfo, view, pfoClusters);
            if (!pfoClusters.empty())
            {
                const Cluster *pCluster{pfoClusters.front()};
                pfoToViewClusterMap[pPfo][view] = pCluster;
                clusterToPfoMap[pCluster] = pPfo;
                if (view != TPC_3D)
                {
                    auto result = this->FitAndOrderCluster(pCluster, clusterToHitMap[pCluster]);
                    if (result)
                        clusterToFitMap.emplace(pCluster, *result);
                }
                else
                    LArClusterHelper::GetAllHits(pCluster, clusterToHitMap[pCluster]);
                for (const CaloHit *const pCaloHit : clusterToHitMap[pCluster])
                {
                    hitToClusterMap[pCaloHit] = pCluster;
                    viewToHitsMap[view].emplace_back(pCaloHit);
                }
            }
        }
    }

    // Loop over the PFOs and project 3D hits from two views into the third view to look for unassociated hits
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    for (const Pfo *const pPfo : *pPfoList)
    {
        const Cluster *const pClusterU{pfoToViewClusterMap[pPfo][TPC_VIEW_U]}, *const pClusterV{pfoToViewClusterMap[pPfo][TPC_VIEW_V]},
            *const pClusterW{pfoToViewClusterMap[pPfo][TPC_VIEW_W]};
        const CaloHitList &hitsU{clusterToHitMap[pClusterU]}, &hitsV{clusterToHitMap[pClusterV]}, &hitsW{clusterToHitMap[pClusterW]};

        CaloHitSet unmatchedHitsU, unmatchedHitsV, unmatchedHitsW;
        this->FindUnmatchedHits(hitsU, hitsV, hitsW, unmatchedHitsV, unmatchedHitsW);
        this->FindUnmatchedHits(hitsV, hitsU, hitsW, unmatchedHitsU, unmatchedHitsW);
        this->FindUnmatchedHits(hitsW, hitsU, hitsV, unmatchedHitsU, unmatchedHitsV);

        std::cout << "PFO has " << unmatchedHitsU.size() << " unmatched U hits, "
                  << unmatchedHitsV.size() << " unmatched V hits and " << unmatchedHitsW.size() << " unmatched W hits." << std::endl;

        // Construct 3D hits from the combinatorics of each view pair and project into the third view
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsU, "U", RED));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsV, "V", GREEN));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsW, "W", BLUE));
        for (const CaloHit *const pCaloHit : unmatchedHitsU)
        {
            const CartesianVector position{pCaloHit->GetPositionVector()};
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "Unmatched U", RED, 2.f));
        }
        // HERE: We can use the longitudinal position along the fit relative to the front and back hits to understand if our
        // hits are along the track, before, or after it.
        // If the hits are before or after, we can probably trivial agree to the inclusion/transfer of the hits into this PFO.
        // If the hits are long the track, we can look up the layer, use some precompued layer to hit mappings to find the
        // hits "either side" of the unmatched hit and see if it forms a blocking path - if it does, include it, otherwise don't
        {
            const CartesianVector front{hitsV.front()->GetPositionVector()};
            const CartesianVector back{hitsV.back()->GetPositionVector()};
            float rL{0.f}, rT{0.f};
            if (clusterToFitMap.find(pClusterV) != clusterToFitMap.end())
                clusterToFitMap.at(pClusterV).GetLocalPosition(front, rL, rT);
            std::cout << "Front V hit has local coordinates rL = " << rL << " and rT = " << rT << " in the V cluster frame." << std::endl;
            if (clusterToFitMap.find(pClusterV) != clusterToFitMap.end())
                clusterToFitMap.at(pClusterV).GetLocalPosition(back, rL, rT);
            std::cout << " Back V hit has local coordinates rL = " << rL << " and rT = " << rT << " in the V cluster frame." << std::endl;
        }
        for (const CaloHit *const pCaloHit : unmatchedHitsV)
        {
            const CartesianVector position{pCaloHit->GetPositionVector()};
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "Unmatched V", GREEN, 2.f));
            if (clusterToFitMap.find(pClusterV) == clusterToFitMap.end())
                continue;
            float rL{0.f}, rT{0.f};
            clusterToFitMap.at(pClusterV).GetLocalPosition(position, rL, rT);
            std::cout << "Unmatched U hit has local coordinates rL = " << rL << " and rT = " << rT << " in the V cluster frame." << std::endl;
        }
        for (const CaloHit *const pCaloHit : unmatchedHitsW)
        {
            const CartesianVector position{pCaloHit->GetPositionVector()};
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "Unmatched W", BLUE, 2.f));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::optional<TwoDSlidingFitResult> TrackRecoveryAlgorithm::FitAndOrderCluster(const Cluster *const pCluster, CaloHitList &orderedHits) const
{
    if (pCluster->GetNCaloHits() < 3)
    {
        LArClusterHelper::GetAllHits(pCluster, orderedHits);
        return std::nullopt;
    }
    const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
    TwoDSlidingFitResult sfr(pCluster, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
    LArClusterHelper::OrderHitsAlongTrajectory(pCluster, sfr, orderedHits);

    return sfr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::FindUnmatchedHits(const CaloHitList &hitsA, const CaloHitList &hitsB, const CaloHitList &hitsC, CaloHitSet &unmatchedHitsB,
    CaloHitSet &unmatchedHitsC) const
{
    if (hitsA.empty() || hitsB.empty() || hitsC.empty())
        return;
    const LArPlaneContextObject *pPlaneContextObject{dynamic_cast<const LArPlaneContextObject*>(PandoraContentApi::GetEventContextObject(
        *this, "PlaneContext"))};
    const HitType viewB{hitsB.front()->GetHitType()}, viewC{hitsC.front()->GetHitType()};
    for (const CaloHit *const pCaloHitA : hitsA)
    {
        const LArPlaneContextObject::HitTriplet *pTriplet{pPlaneContextObject->GetHitTriplet(pCaloHitA)};
        if (!pTriplet)
            continue;
        const CaloHit *const pHitB{viewB == TPC_VIEW_W ? pTriplet->m_wHit : viewB == TPC_VIEW_V ? pTriplet->m_vHit : pTriplet->m_uHit};
        const CaloHit *const pHitC{viewC == TPC_VIEW_W ? pTriplet->m_wHit : viewC == TPC_VIEW_V ? pTriplet->m_vHit : pTriplet->m_uHit};
        auto bIterator{std::find(hitsB.begin(), hitsB.end(), pHitB)};
        auto cIterator{std::find(hitsC.begin(), hitsC.end(), pHitC)};
        // Only flag hits as unmatched if the other hits in the triplet are present in the PFO
        if (bIterator == hitsB.end() && cIterator != hitsC.end())
            unmatchedHitsB.insert(pHitB);
        if (cIterator == hitsC.end() && bIterator != hitsB.end())
            unmatchedHitsC.insert(pHitC);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
