/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include <tuple>

#include <Eigen/Dense>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

std::atomic<int> KalmanClusterCreationAlgorithm::KalmanFit::m_counter(0);

KalmanClusterCreationAlgorithm::KalmanClusterCreationAlgorithm() :
    m_kalmanDelta(1.f),
    m_kalmanProcessVarCoeff(1.f),
    m_kalmanMeasurementVarCoeff(1.f),
    m_minCosTheta(0.94f),
    m_xTol(0.25f),
    m_zTol(0.125f),
    m_kalmanMinHits(4)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanClusterCreationAlgorithm::~KalmanClusterCreationAlgorithm()
{
    this->Reset();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::Run()
{
    this->Reset();
    // Collect hits into layer-ordered hit lists
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;
    HitType view{pCaloHitList->front()->GetHitType()};
    m_viewHitsMap[view] = OrderedCaloHitList();
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FilterCaloHits(pCaloHitList, m_viewHitsMap[view]));

    float min{std::numeric_limits<float>::max()}, max{std::numeric_limits<float>::lowest()};
    this->GetSpanX(min, max);
    CaloHitList caloHitList;
    m_viewHitsMap[view].FillCaloHitList(caloHitList);
    if (caloHitList.empty())
        return STATUS_CODE_SUCCESS;
    m_slicedCaloHits[view] = new LArSlicedCaloHitList(caloHitList, min, max);
    this->IdentifyCandidateClusters(view);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::IdentifyCandidateClusters(const HitType view)
{
    std::map<HitType, IDKalmanFitMap> viewKalmanFitsMap;
    std::map<HitType, HitKalmanFitMap> viewHitKalmanFitMap;
    viewKalmanFitsMap[view] = IDKalmanFitMap();
    viewHitKalmanFitMap[view] = HitKalmanFitMap();
    IDKalmanFitMap &kalmanFits{viewKalmanFitsMap[view]};
    HitKalmanFitMap &hitKalmanFitMap{viewHitKalmanFitMap[view]};

    for (const auto &[bin, caloHits0] : m_slicedCaloHits[view]->GetSliceHitMap())
    {
        this->MakeClusterSeeds(caloHits0, kalmanFits, hitKalmanFitMap);
        this->BuildClusters(caloHits0, kalmanFits, hitKalmanFitMap);
        this->RemoveDuplicateKalmanFits(kalmanFits, hitKalmanFitMap);
        this->AllocateAmbiguousHits(kalmanFits, hitKalmanFitMap);
        this->ConsolidateInterleavedClusters(caloHits0, kalmanFits, hitKalmanFitMap);
    }

    /*PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    for (auto &[id, kalmanFit] : kalmanFits)
    {
        const CaloHitList caloHits(kalmanFit.m_caloHits.begin(), kalmanFit.m_caloHits.end());
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "Kalman", AUTOITER));
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/

    const ClusterList *pClusterList{nullptr};
    std::string temporaryListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, temporaryListName));
    for (auto &[id, kalmanFit] : kalmanFits)
    {
        const Cluster *pCluster{nullptr};
        for(const CaloHit *const pCaloHit : kalmanFit.m_caloHits)
        {
            if (pCluster)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit));
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
            }
        }
    }
    const std::string viewStr{view == TPC_VIEW_U ? "U" : view == TPC_VIEW_V ? "V" : "W"};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListPrefix + viewStr));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::MakeClusterSeeds(const CaloHitVector &sliceCaloHits, IDKalmanFitMap &kalmanFits, HitKalmanFitMap &hitKalmanFitMap)
{
    if (sliceCaloHits.empty())
        return;
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view{sliceCaloHits.front()->GetHitType()};
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};

    for (auto iter1 = sliceCaloHits.begin(); iter1 != sliceCaloHits.end(); ++iter1)
    {
        const CaloHit *const pSeedHit{*iter1};
        const CartesianVector &seed{pSeedHit->GetPositionVector()};
        Eigen::VectorXd init(2);
        init << seed.GetX(), seed.GetZ();
        for (auto iter2 = std::next(iter1); iter2 != sliceCaloHits.end(); ++iter2)
        {
            const CaloHit *const pCaloHit{*iter2};

            const CaloHitList caloHits(sliceCaloHits.begin(), sliceCaloHits.end());
            const CartesianVector &other{pCaloHit->GetPositionVector()};
            if (this->Proximate(pSeedHit, pCaloHit))
            {
                if (this->SkipsOverHit(sliceCaloHits, pSeedHit, pCaloHit))
                    continue;
                KalmanFit kalmanFit{pSeedHit, KalmanFilter2D(m_kalmanDelta, processVariance, measurementVariance, init), CaloHitVector()};
                
                kalmanFit.m_kalmanFilter.Predict();
                Eigen::VectorXd measurement(2);
                measurement << other.GetX(), other.GetZ();
                kalmanFit.m_kalmanFilter.Update(measurement);
                kalmanFit.InsertHit(pSeedHit);
                hitKalmanFitMap[pSeedHit].insert(kalmanFit.m_id);
                kalmanFit.InsertHit(pCaloHit);
                hitKalmanFitMap[pCaloHit].insert(kalmanFit.m_id);
                kalmanFits.emplace(kalmanFit.m_id, kalmanFit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::BuildClusters(const CaloHitVector &sliceCaloHits, IDKalmanFitMap &kalmanFits, HitKalmanFitMap &hitKalmanFitMap)
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    for (auto &[id, kalmanFit] : kalmanFits)
    {
        bool added{true};
        while (added)
        {
            added = false;
            kalmanFit.m_kalmanFilter.Predict();
            const Eigen::VectorXd &state{kalmanFit.m_kalmanFilter.GetTemporaryState()};
            Eigen::VectorXd previousMeasurement(2);
            previousMeasurement << kalmanFit.m_pLastHit->GetPositionVector().GetX(), kalmanFit.m_pLastHit->GetPositionVector().GetZ();

            const CaloHit *pBestHit{nullptr};
            Eigen::VectorXd bestMeasurement(2);
            double bestDistanceSquared{std::numeric_limits<double>::max()};
            for (const CaloHit *const pCaloHit : sliceCaloHits)
            {
                if (std::find(kalmanFit.m_caloHits.begin(), kalmanFit.m_caloHits.end(), pCaloHit) != kalmanFit.m_caloHits.end())
                    continue;
                const CartesianVector &other{pCaloHit->GetPositionVector()};
                if ((kalmanFit.m_caloHits.size() < 4 && this->Proximate(kalmanFit.m_pLastHit, pCaloHit)) || this->Contains(pCaloHit, state, m_xTol, m_zTol))
                {
                    Eigen::VectorXd measurement(2);
                    measurement << other.GetX(), other.GetZ();
                    // Check that the hit is (very approximately) in the direction of the Kalman filter
                    const Eigen::VectorXd &thisDirection{(measurement - previousMeasurement).normalized()};
                    const double cosTheta{thisDirection.dot(kalmanFit.m_kalmanFilter.GetDirection())};

                    if (false && pCaloHit->GetHitType() == TPC_VIEW_V)
                    {
                        const CaloHitList caloHits(kalmanFit.m_caloHits.begin(), kalmanFit.m_caloHits.end());
                        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "Input", GRAY));
                        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &other, "Measurement", BLUE, 2));
                        const CartesianVector predicted{CartesianVector(state[0], 0.f, state[1])};
                        const CartesianVector lastHit{kalmanFit.m_pLastHit->GetPositionVector()};
                        const CartesianVector direction{lastHit + CartesianVector(kalmanFit.m_kalmanFilter.GetDirection()[0], 0.f, kalmanFit.m_kalmanFilter.GetDirection()[1]) * 5};
                        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &predicted, "Prediction " + std::to_string((state - measurement).norm()), RED, 2));
                        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &lastHit, &direction, "Direction " + std::to_string(cosTheta), BLACK, 1, 1));
                        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                    }

                    if (cosTheta < m_minCosTheta)
                        continue;
                    const double distanceSquared{(state - measurement).squaredNorm()};
                    if (distanceSquared < bestDistanceSquared)
                    {
                        bestDistanceSquared = distanceSquared;
                        bestMeasurement = measurement;
                        pBestHit = pCaloHit;
                    }
                }
            }
            if (pBestHit)
            {
                added = true;
                kalmanFit.m_kalmanFilter.Update(bestMeasurement);
                kalmanFit.InsertHit(pBestHit);
                hitKalmanFitMap[pBestHit].insert(kalmanFit.m_id);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::RemoveDuplicateKalmanFits(IDKalmanFitMap &kalmanFits, HitKalmanFitMap &hitKalmanFitMap)
{
    for (auto iter1 = kalmanFits.begin(); iter1 != kalmanFits.end(); ++iter1)
    {
        KalmanFit &kalmanFit1{iter1->second};
        for (auto iter2 = std::next(iter1); iter2 != kalmanFits.end();)
        {
            KalmanFit &kalmanFit2{iter2->second};
            bool subset{true};
            for (const CaloHit *const pCaloHit : kalmanFit2.m_caloHits)
            {
                if (hitKalmanFitMap[pCaloHit].find(kalmanFit1.m_id) == hitKalmanFitMap[pCaloHit].end())
                {
                    subset = false;
                    break;
                }
            }
            if (subset)
            {
                for (const CaloHit *const pCaloHit : kalmanFit2.m_caloHits)
                {
                    auto found{hitKalmanFitMap[pCaloHit].find(kalmanFit2.m_id)};
                    if (found != hitKalmanFitMap[pCaloHit].end())
                        hitKalmanFitMap[pCaloHit].erase(found);
                }
                iter2 = kalmanFits.erase(iter2);
            }
            else
            {
                ++iter2;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::AllocateAmbiguousHits(IDKalmanFitMap &kalmanFits, HitKalmanFitMap &hitKalmanFitMap)
{
    HitKalmanFitMap removals;
    for (const auto &[pCaloHit, fitIDs] : hitKalmanFitMap)
    {
        if (fitIDs.size() == 1)
            continue;
        std::map<int, double> idErrorMap;
        for (const int id : fitIDs)
        {
            if (kalmanFits.find(id) != kalmanFits.end())
            {
                const CaloHitVector &caloHits{kalmanFits.at(id).m_caloHits};
                auto targetHit{std::find(caloHits.begin(), caloHits.end(), pCaloHit)};
                if (targetHit != caloHits.end())
                {
                    Eigen::VectorXd init(2);
                    init << caloHits.front()->GetPositionVector().GetX(), caloHits.front()->GetPositionVector().GetZ();
                    Eigen::VectorXd target(2);
                    target << (*targetHit)->GetPositionVector().GetX(), (*targetHit)->GetPositionVector().GetZ();
                    const auto dForward{std::distance(caloHits.begin(), targetHit)};
                    const auto dBackward{static_cast<long>(caloHits.size()) - dForward - 1};
                    // Look at the quality of the different Kalman fits for ambiguous hits
                    if (dForward >= 2 && dBackward >= 2)
                    {
                        // Either enough hits to check fit quality from both directions, or it's a small cluster so we do the best
                        // we can with limited information
                        Eigen::VectorXd last(2);
                        last << caloHits.back()->GetPositionVector().GetX(), caloHits.back()->GetPositionVector().GetZ();
                        idErrorMap[id] = std::min(this->GetForwardError(caloHits, init, target, dForward),
                            this->GetBackwardError(caloHits, last, target, dBackward));
                    }
                    else if (dForward >= 2)
                    {
                        // Fit more reliable in forward direction, use that
                        idErrorMap[id] = this->GetForwardError(caloHits, init, target, dForward);
                    }
                    else if (dBackward >= 2)
                    {
                        // Fit more reliable in backward direction, use that
                        Eigen::VectorXd last(2);
                        last << caloHits.back()->GetPositionVector().GetX(), caloHits.back()->GetPositionVector().GetZ();
                        idErrorMap[id] = this->GetBackwardError(caloHits, last, target, dBackward);
                    }
                    else
                    {
                        // 2 hit cluster - apply penalty
                        idErrorMap[id] = 1000;
                    }
                }
            }
        }
        // Find the best match, but dont update yet, as it could have side-effects on the subsequent Kalman filters due to "missing" hits
        double bestError{std::numeric_limits<double>::max()};
        int bestId{-1};
        for (const auto &[id, error] : idErrorMap)
        {
            if (error < bestError)
            {
                bestError = error;
                bestId = id;
            }
        }

        // Update the hit to fit map to only retain the best cluster and retain the other ids for ease of removal later
        if (bestId > -1)
        {
            std::copy_if(hitKalmanFitMap.at(pCaloHit).begin(), hitKalmanFitMap.at(pCaloHit).end(), std::inserter(removals[pCaloHit],
                removals[pCaloHit].begin()), [&bestId](int id){ return id != bestId; });
            hitKalmanFitMap[pCaloHit] = {bestId};
        }
    }
    for (const auto &[id, kalmanFit] : kalmanFits)
    {
        IntVector competingIds;
        for (const CaloHit *const pCaloHit : kalmanFit.m_caloHits)
        {
            if (hitKalmanFitMap.find(pCaloHit) != hitKalmanFitMap.end())
            {
                for (const int i : hitKalmanFitMap.at(pCaloHit))
                    if (id != i)
                        competingIds.emplace_back(i);
            }
        }
        std::sort(competingIds.begin(), competingIds.end());
        for (const int i : competingIds)
        {
            // Only consider reallocation of hits now if this is the larger cluster, or, if they're the same size, if this is the earlier cluster
            if ((kalmanFit.m_caloHits.size() < kalmanFits.at(i).m_caloHits.size()) ||
                (id > i && kalmanFit.m_caloHits.size() == kalmanFits.at(i).m_caloHits.size()))
                continue;

            CaloHitVector intersection;
            std::set_intersection(kalmanFit.m_caloHits.begin(), kalmanFit.m_caloHits.end(), kalmanFits.at(i).m_caloHits.begin(),
                kalmanFits.at(i).m_caloHits.end(), std::back_inserter(intersection));

            if (intersection.size() > 1)
            {
                for (const CaloHit *const pCaloHit : intersection)
                {
                    if (hitKalmanFitMap.at(pCaloHit).find(i) != hitKalmanFitMap.at(pCaloHit).end())
                    {
                        hitKalmanFitMap[pCaloHit] = {id};
                        removals[pCaloHit].insert(i);
                        auto found{removals[pCaloHit].find(id)};
                        if (found != removals[pCaloHit].end())
                            removals[pCaloHit].erase(found);
                    }
                }
            }
        }
    }

    // Now we're done with finding the best cluster for each hit, remove the respective hits from other clusters
    for (const auto &[pCaloHit, ids] : removals)
    {
        for (const int id : ids)
        {
            if (kalmanFits.find(id) != kalmanFits.end())
            {
                CaloHitVector &caloHits{kalmanFits.at(id).m_caloHits};
                auto found{std::find(caloHits.begin(), caloHits.end(), pCaloHit)};
                if (found != caloHits.end())
                {
                    caloHits.erase(found);
                    if (hitKalmanFitMap.find(pCaloHit) != hitKalmanFitMap.end())
                    {
                        auto found2{std::find(hitKalmanFitMap.at(pCaloHit).begin(), hitKalmanFitMap.at(pCaloHit).end(), id)};
                        if (found2 != hitKalmanFitMap.at(pCaloHit).end())
                            hitKalmanFitMap.at(pCaloHit).erase(found2);
                    }
                }
            }
        }
    }

    // Remove any one hit clusters
    for (auto iter = kalmanFits.begin(); iter != kalmanFits.end(); )
    {
        const int id{iter->second.m_id};
        if (iter->second.m_caloHits.size() <= 1)
        {
            for (const CaloHit *const pCaloHit : iter->second.m_caloHits)
            {
                if (hitKalmanFitMap.find(pCaloHit) != hitKalmanFitMap.end())
                {
                    auto found{std::find(hitKalmanFitMap.at(pCaloHit).begin(), hitKalmanFitMap.at(pCaloHit).end(), id)};
                    if (found != hitKalmanFitMap.at(pCaloHit).end())
                        hitKalmanFitMap.at(pCaloHit).erase(found);
                }
            }
            iter = kalmanFits.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::ConsolidateInterleavedClusters(const CaloHitVector &caloHits, IDKalmanFitMap &kalmanFits, HitKalmanFitMap &hitKalmanFitMap)
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    KalmanFitIDSet sliceIDs;
    // Identify the fits that are part of this slice
    for (const CaloHit *const pCaloHit : caloHits)
    {
        if (hitKalmanFitMap.find(pCaloHit) != hitKalmanFitMap.end())
        {
            for (const int id : hitKalmanFitMap.at(pCaloHit))
                sliceIDs.insert(id);
        }
    }
    if (sliceIDs.size() < 2)
        return;
    std::vector<std::pair<int, int>> sliceIDHitPairs;
    for (const int id : sliceIDs)
        sliceIDHitPairs.emplace_back(std::make_pair(kalmanFits.at(id).m_caloHits.size(), id));
    if (sliceIDHitPairs.size() < 2)
        return;
    std::sort(sliceIDHitPairs.begin(), sliceIDHitPairs.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.first > b.first; });

    std::map<int, CaloHitVector> hitVectors;
    for (auto &[hitCount, id] : sliceIDHitPairs)
    {
        const CaloHitVector &fitCaloHits{kalmanFits.at(id).m_caloHits};
        for (const CaloHit *const pCaloHit : fitCaloHits)
        {
            if (std::find(caloHits.begin(), caloHits.end(), pCaloHit) != caloHits.end())
                hitVectors[id].emplace_back(pCaloHit);
        }
    }

    for (auto pairIter1 = sliceIDHitPairs.begin(); pairIter1 != sliceIDHitPairs.end(); ++pairIter1)
    {
        const CaloHitVector &hits1{hitVectors[pairIter1->second]};
        std::vector<std::tuple<int, const CaloHit *, const CaloHit *>> collectedHits;
        // Get consecutive hits in this cluster and look for intervening hits in the other
        for (auto hIter = hits1.begin(); hIter != hits1.end(); ++hIter)
        {
            const CaloHit *const pCaloHit1{*hIter};
            if (std::next(hIter) == hits1.end())
                break;
            const CaloHit *const pCaloHit2{*std::next(hIter)};
            for (auto pairIter2 = std::next(pairIter1); pairIter2 != sliceIDHitPairs.end(); ++pairIter2)
            {
                const CaloHitVector &hits2{hitVectors[pairIter2->second]};
                for (const CaloHit *const pCheckedHit : hits2)
                {
                    CaloHitVector singleHit{pCheckedHit};
                    if (this->SkipsOverHit(singleHit, pCaloHit1, pCaloHit2))
                        collectedHits.emplace_back(std::make_tuple(pairIter2->second, pCaloHit2, pCheckedHit));
                }
            }
        }
        if (collectedHits.empty())
            continue;

/*        const CaloHitList caloHits1(hits1.begin(), hits1.end());
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits1, "Main", BLUE));
        CaloHitList caloHits2;
        for (const auto &[id, pCaloHit1, pCaloHit2] : collectedHits)
            caloHits2.emplace_back(pCaloHit2);
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits2, "Collected", RED));*/

        // There is something wrong here, the hits are not being moved correctly
        for (const auto &[id, pCaloHit1, pCaloHit2] : collectedHits)
        {
            CaloHitVector &otherHits{kalmanFits.at(id).m_caloHits};
            auto oldPos{std::find(otherHits.begin(), otherHits.end(), pCaloHit2)};
            // There are edge cases where a hit can be between more than one pair of hits, need to skip these to avoid double removal
            if (oldPos == otherHits.end())
                continue;
            otherHits.erase(oldPos);
            oldPos = std::find(hitVectors.at(id).begin(), hitVectors.at(id).end(), pCaloHit2);
            hitVectors.at(id).erase(oldPos);
            CaloHitVector &currentHits{kalmanFits.at(pairIter1->second).m_caloHits};
            auto newPos{std::find(currentHits.begin(), currentHits.end(), pCaloHit1)};
            currentHits.insert(newPos, pCaloHit2);
            newPos = std::find(hitVectors.at(pairIter1->second).begin(), hitVectors.at(pairIter1->second).end(), pCaloHit1);
            hitVectors.at(pairIter1->second).insert(newPos, pCaloHit2);
            hitKalmanFitMap[pCaloHit2] = {pairIter1->second};
        }
        /*
        const CaloHitList caloHits3(kalmanFits.at(vIter1->first).m_caloHits.begin(), kalmanFits.at(vIter1->first).m_caloHits.end());
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits3, "Updated", BLACK));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
    }
    // Remove any one hit clusters
    for (auto iter = kalmanFits.begin(); iter != kalmanFits.end(); )
    {
        const int id{iter->second.m_id};
        if (iter->second.m_caloHits.size() <= 1)
        {
            for (const CaloHit *const pCaloHit : iter->second.m_caloHits)
            {
                if (hitKalmanFitMap.find(pCaloHit) != hitKalmanFitMap.end())
                {
                    auto found{std::find(hitKalmanFitMap.at(pCaloHit).begin(), hitKalmanFitMap.at(pCaloHit).end(), id)};
                    if (found != hitKalmanFitMap.at(pCaloHit).end())
                        hitKalmanFitMap.at(pCaloHit).erase(found);
                }
            }
            iter = kalmanFits.erase(iter);
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

double KalmanClusterCreationAlgorithm::GetForwardError(const CaloHitVector &caloHits, const Eigen::VectorXd &init, const Eigen::VectorXd &target,
    const long steps)
{
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    // caloHits cannot be empty, as empty Kalman fits are deleted
    const HitType view{caloHits.front()->GetHitType()};
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};

    KalmanFilter2D kalman{m_kalmanDelta, processVariance, measurementVariance, init};
    auto iter{std::next(caloHits.begin())};
    for (long i = 1; i <= steps; ++i, ++iter)
    {
        kalman.Predict();
        const CartesianVector pos{(*iter)->GetPositionVector()};
        Eigen::VectorXd measurement(2);
        measurement << pos.GetX(), pos.GetZ();
        kalman.Update(measurement);
    }
    return (target - kalman.GetPosition()).squaredNorm();
}

//------------------------------------------------------------------------------------------------------------------------------------------

double KalmanClusterCreationAlgorithm::GetBackwardError(const CaloHitVector &caloHits, const Eigen::VectorXd &init, const Eigen::VectorXd &target,
    const long steps)
{
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    // caloHits cannot be empty, as empty Kalman fits are deleted
    const HitType view{caloHits.front()->GetHitType()};
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};

    KalmanFilter2D kalman{m_kalmanDelta, processVariance, measurementVariance, init};
    auto iter{std::next(caloHits.rbegin())};
    for (long i = 1; i <= steps; ++i, ++iter)
    {
        kalman.Predict();
        const CartesianVector pos{(*iter)->GetPositionVector()};
        Eigen::VectorXd measurement(2);
        measurement << pos.GetX(), pos.GetZ();
        kalman.Update(measurement);
    }
    return (target - kalman.GetPosition()).squaredNorm();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::Make3DHitPermutations(const CaloHitVector &caloHits0, const CaloHitVector &caloHits1, const CaloHitVector &caloHits2, CandidateCluster::HitTripletVector &triplets) const
{
    int i{0};
    // Construct 3 hit tuples
    for (const CaloHit *pCaloHit0 : caloHits0)
    {
        for (const CaloHit *pCaloHit1 : caloHits1)
        {
            for (const CaloHit *pCaloHit2 : caloHits2)
            {
                std::get<0>(triplets[i]) = pCaloHit0;
                std::get<1>(triplets[i]) = pCaloHit1;
                std::get<2>(triplets[i]) = pCaloHit2;
                ++i;
            }
        }
    }
    // Construct 2 hit tuples
    for (const CaloHit *pCaloHit0 : caloHits0)
    {
        for (const CaloHit *pCaloHit1 : caloHits1)
        {
            std::get<0>(triplets[i]) = pCaloHit0;
            std::get<1>(triplets[i]) = pCaloHit1;
            ++i;
        }
    }
    for (const CaloHit *pCaloHit0 : caloHits0)
    {
        for (const CaloHit *pCaloHit2 : caloHits2)
        {
            std::get<0>(triplets[i]) = pCaloHit0;
            std::get<1>(triplets[i]) = pCaloHit2;
            ++i;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::FilterCaloHits(const CaloHitList *const pCaloHitList, OrderedCaloHitList &selectedCaloHitList) const
{
    CaloHitList availableHitList;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (PandoraContentApi::IsAvailable(*this, pCaloHit))
            availableHitList.push_back(pCaloHit);
    }

    if (availableHitList.empty())
        return STATUS_CODE_SUCCESS;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, selectedCaloHitList.Add(availableHitList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::GetSpanX(float &min, float &max) const
{
    min = std::numeric_limits<float>::max();
    max = std::numeric_limits<float>::lowest();

    for (const auto &[view, orderedHitList] : m_viewHitsMap)
    {
        for (OrderedCaloHitList::const_iterator iter = orderedHitList.begin(), iterEnd = orderedHitList.end(); iter != iterEnd; ++iter)
        {
            const CaloHitList &caloHitList{*iter->second};
            for (const CaloHit *const pCaloHit : caloHitList)
            {
                const float xLo{pCaloHit->GetPositionVector().GetX() - 0.5f * pCaloHit->GetCellSize1()};
                const float xHi{pCaloHit->GetPositionVector().GetX() + 0.5f * pCaloHit->GetCellSize1()};
                if (xLo < min)
                    min = xLo;
                if (xHi > max)
                    max = xHi;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::Reset()
{
    StatusCode status{Algorithm::Reset()};
    for (auto &[bin, slicedCaloHitList] : m_slicedCaloHits)
    {
        delete slicedCaloHitList;
    }
    m_slicedCaloHits.clear();
    m_viewHitsMap.clear();

    return status;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::GetSlices(const LArSlicedCaloHitList::SliceHitMap &sliceHitMap, const size_t bin, CaloHitVector &caloHits_m1,
    CaloHitVector &caloHits_p1) const
{
    if (sliceHitMap.find(bin - 1) != sliceHitMap.end())
    {
        const CaloHitVector &caloHits{sliceHitMap.at(bin - 1)};
        caloHits_m1.insert(caloHits_m1.end(), caloHits.begin(), caloHits.end());
    }
    if (sliceHitMap.find(bin + 1) != sliceHitMap.end())
    {
        const CaloHitVector &caloHits{sliceHitMap.at(bin + 1)};
        caloHits_p1.insert(caloHits_p1.end(), caloHits.begin(), caloHits.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::GetSlices(const LArSlicedCaloHitList::SliceHitMap &sliceHitMap, const size_t bin, CaloHitVector &caloHits) const
{
    for (int i = -1; i <= 1; ++i)
    {
        if (sliceHitMap.find(bin + i) != sliceHitMap.end())
        {
            const CaloHitVector &sliceCaloHits{sliceHitMap.at(bin + i)};
            caloHits.insert(caloHits.end(), sliceCaloHits.begin(), sliceCaloHits.end());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanDelta", m_kalmanDelta));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanProcessVarCoeff", m_kalmanProcessVarCoeff));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanMeasurementVarCoeff", m_kalmanMeasurementVarCoeff));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "m_minCosTheta", m_minCosTheta));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "m_xTol", m_xTol));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "m_zTol", m_zTol));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanMinHits", m_kalmanMinHits));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListPrefix", m_clusterListPrefix));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KalmanClusterCreationAlgorithm::Proximate(const CaloHit *const pCaloHit1, const CaloHit *const pCaloHit2, const float proximity) const
{
    const CartesianVector &position1{pCaloHit1->GetPositionVector()};
    const CartesianVector &position2{pCaloHit2->GetPositionVector()};

    const float distanceSquared{(position1 - position2).GetMagnitudeSquared()};
    if (distanceSquared < proximity * proximity)
        return true;

    if (std::fabs(position1.GetZ() - position2.GetZ()) >= proximity)
        return false;

    const float width1{0.5f * pCaloHit1->GetCellSize1()};
    const float width2{0.5f * pCaloHit2->GetCellSize1()};
    const float xlo1{position1.GetX() - width1 - proximity}, xhi1{position1.GetX() + width1 + proximity};
    const float xlo2{position2.GetX() - width2 - proximity}, xhi2{position2.GetX() + width2 + proximity};
    return ((position1.GetX() >= xlo2 && position1.GetX() <= xhi2) || (position2.GetX() >= xlo1 && position2.GetX() <= xhi1));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KalmanClusterCreationAlgorithm::Contains(const CaloHit *const pCaloHit, const Eigen::VectorXd &position, const float xTol, const float zTol) const
{
    const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
    const float width{0.5f * pCaloHit->GetCellSize1() + xTol};
    const float xlo{hitPosition.GetX() - width}, xhi{hitPosition.GetX() + width};
    const float height{0.5f * pCaloHit->GetCellSize0() + zTol};
    const float zlo{hitPosition.GetZ() - height}, zhi{hitPosition.GetZ() + height};
    return (position(0) >= xlo && position(0) <= xhi && position(1) >= zlo && position(1) <= zhi);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KalmanClusterCreationAlgorithm::SkipsOverHit(const CaloHitVector &caloHits, const CaloHit *const pCaloHit1, const CaloHit *const pCaloHit2) const
{
    // For each hit in the calo hit list, check if the line between the two hits passes through the hit
    const CartesianVector &pos1{pCaloHit1->GetPositionVector()};
    const CartesianVector &pos2{pCaloHit2->GetPositionVector()};
    for (const CaloHit *const pCaloHit : caloHits)
    {
        if (pCaloHit == pCaloHit1 || pCaloHit == pCaloHit2)
            continue;

        const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
        const double xmin{hitPosition.GetX() - 0.5 * pCaloHit->GetCellSize1()}, xmax{hitPosition.GetX() + 0.5 * pCaloHit->GetCellSize1()};
        const double zmin{hitPosition.GetZ() - 0.5 * pCaloHit->GetCellSize0()}, zmax{hitPosition.GetZ() + 0.5 * pCaloHit->GetCellSize0()};
        const double hitSize{(xmax - xmin) * (zmax - zmin)};

        double entry{0.}, exit{1.};

        auto check_axis = [&](double p1, double p2, double minB, double maxB)
        {
            double t1{(minB - p1) / (p2 - p1)};
            double t2{(maxB - p1) / (p2 - p1)};

            if (t1 > t2)
                std::swap(t1, t2);

            entry = std::max(entry, t1);
            exit = std::min(exit, t2);

            return entry <= exit;
        };

        if (check_axis(pos1.GetX(), pos2.GetX(), xmin, xmax) && check_axis(pos1.GetZ(), pos2.GetZ(), zmin, zmax))
        {
            // Check if the intervening hit and the first extremal hit have a large overlap
            const double xmin1{pos1.GetX() - 0.5 * pCaloHit1->GetCellSize1()}, xmax1{pos1.GetX() + 0.5 * pCaloHit1->GetCellSize1()};
            const double zmin1{pos1.GetZ() - 0.5 * pCaloHit1->GetCellSize0()}, zmax1{pos1.GetZ() + 0.5 * pCaloHit1->GetCellSize0()};
            const double hitSize1{(xmax1 - xmin1) * (zmax1 - zmin1)};
            const double x_left1{std::max(xmin, xmin1)}, x_right1{std::min(xmax, xmax1)};
            const double z_bottom1{std::max(zmin, zmin1)}, z_top1{std::min(zmax, zmax1)};
            double overlap1{0};
            if (x_left1 < x_right1 && z_bottom1 < z_top1)
                overlap1 = (x_right1 - x_left1) * (z_top1 - z_bottom1);
            if (overlap1 > 0.5 * hitSize1 || overlap1 > 0.5 * hitSize)
                return false;

            // Check if the intervening hit and the second extremal hit have a large overlap
            const double xmin2{pos2.GetX() - 0.5 * pCaloHit2->GetCellSize1()}, xmax2{pos2.GetX() + 0.5 * pCaloHit2->GetCellSize1()};
            const double zmin2{pos2.GetZ() - 0.5 * pCaloHit2->GetCellSize0()}, zmax2{pos2.GetZ() + 0.5 * pCaloHit2->GetCellSize0()};
            const double hitSize2{(xmax1 - xmin1) * (zmax1 - zmin1)};
            const double x_left2{std::max(xmin, xmin2)}, x_right2{std::min(xmax, xmax2)};
            const double z_bottom2{std::max(zmin, zmin2)}, z_top2{std::min(zmax, zmax2)};
            double overlap2{0};
            if (x_left2 < x_right2 && z_bottom2 < z_top2)
                overlap2 = (x_right2 - x_left2) * (z_top2 - z_bottom2);
            if (overlap2 > 0.5 * hitSize2 || overlap2 > 0.5 * hitSize)
                return false;

            // We have an intervening hit without significant overlap with either extremal hit
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::KalmanFit::InsertHit(const CaloHit *const pCaloHit)
{
    m_caloHits.emplace_back(pCaloHit);
    m_pLastHit = pCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet(const CaloHit *const pHit1, const CaloHit *const pHit2, const CaloHit *const pHit3)
{
    const CaloHit *pHitU{nullptr}, *pHitV{nullptr}, *pHitW{nullptr};
    if (pHit1)
    {
        switch (pHit1->GetHitType())
        {
            case TPC_VIEW_U:
                pHitU = pHit1;
                break;
            case TPC_VIEW_V:
                pHitU = pHit1;
                break;
            case TPC_VIEW_W:
                pHitU = pHit1;
                break;
            default:
                std::cout << "Error: KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet - Unexpetced hit type: " << pHit1->GetHitType() << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }
    if (pHit2)
    {
        switch (pHit2->GetHitType())
        {
            case TPC_VIEW_U:
                pHitU = pHit2;
                break;
            case TPC_VIEW_V:
                pHitU = pHit2;
                break;
            case TPC_VIEW_W:
                pHitU = pHit2;
                break;
            default:
                std::cout << "Error: KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet - Unexpetced hit type: " << pHit2->GetHitType() << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }
    if (pHit3)
    {
        switch (pHit3->GetHitType())
        {
            case TPC_VIEW_U:
                pHitU = pHit3;
                break;
            case TPC_VIEW_V:
                pHitU = pHit3;
                break;
            case TPC_VIEW_W:
                pHitU = pHit3;
                break;
            default:
                std::cout << "Error: KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet - Unexpetced hit type: " << pHit3->GetHitType() << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }
    m_hitTriplets.emplace_back(std::make_tuple(pHitU, pHitV, pHitW));
}

} // namespace lar_content
