/**
 *  @file   larpandoracontent/LArMonitoring/TrackOverlapMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TrackOverlapMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

TrackOverlapMonitoringAlgorithm::TrackOverlapMonitoringAlgorithm() :
    m_visualise{true},
    m_writeFile{false},
    m_vertexRadius{10.f},
    m_distance{2.f},
    m_delta{0.5f},
    m_processVarCoeff{1.f},
    m_measurementVarCoeff{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapMonitoringAlgorithm::~TrackOverlapMonitoringAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootFileName, m_rootTreeName, "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::Run()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    }

    this->CreateMaps();
    MCToMCMap overlapCandidates;
    this->FindTrueOverlapCandidates(overlapCandidates);
    std::cout << "Overlap candidates found: " << overlapCandidates.size() << std::endl;
    this->AssessPfos(overlapCandidates);

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::CreateMaps()
{
    const PfoList *pfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pfoList));
    if (!pfoList || pfoList->empty())
        return STATUS_CODE_SUCCESS;

    for (const Pfo *const pPfo : *pfoList)
    {
        CaloHitList allCaloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, allCaloHits);
        for (const CaloHit *const pCaloHit : allCaloHits)
        {
            const MCParticleWeightMap &contributionMap{pCaloHit->GetMCParticleWeightMap()};
            for (const auto &[pMC, weight] : contributionMap)
            {
                if (weight <= 0.f)
                    continue;

                m_pfoToMCMap[pPfo].insert(pMC);
                m_mcToPfoMap[pMC].insert(pPfo);
                m_mcToHitsMap[pMC].insert(pCaloHit);
            }
        }
    }

    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        const int pdg{std::abs(pMC->GetParticleId())};
        if (caloHits.empty() || (pdg == 11|| pdg == 22) || pMC->GetMomentum().GetMagnitude() < 1e-3)
            continue;
        if (m_vertexToMCMap.find(pMC->GetVertex()) == m_vertexToMCMap.end())
            m_vertexToMCMap[pMC->GetVertex()] = MCParticleSet();
        m_vertexToMCMap[pMC->GetVertex()].insert(pMC);
        // Might want to think about endpoints too, but that makes direction comparison harder
        //m_vertexToMCMap[pMC->GetEndPoint()].insert(pMC);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::FindTrueOverlapCandidates(MCToMCMap &overlapCandidates) const
{
    const LArTransformationPlugin *pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (const auto &[vertex, mcParticles] : m_vertexToMCMap)
    {
        if (mcParticles.size() < 2)
            continue;

        // Get all of the hits in each MC particle and find the closest approach for each hit in the other MC particles.
        // Retain hits within a couple of cm of a hit in each of the other MC particles.
        // If a sufficiently large fraction of hits are close within a 10cm radius of the vertex, we have an overlap candidate.
        for (auto iter1 = mcParticles.begin(); iter1 != mcParticles.end(); ++iter1)
        {
            const MCParticle *const pMC1{*iter1};
            CaloHitList uHits1, vHits1, wHits1;
            this->CollectHitsByView(pMC1, uHits1, vHits1, wHits1);
            float x{0}, u{0}, v{0}, w{0};
            LArVertexHelper::GetPositionProjections(pMC1->GetVertex(), pTransform, x, u, v, w);
            Eigen::RowVector2f uVertex(x, u);
            Eigen::MatrixXf uHitMatrix1; // Size will be set by FilterHits
            this->VectorizeAndFilterHits(uHits1, uVertex, m_vertexRadius, uHitMatrix1);
            // Need to do the same filtering for the second batch of hits

            for (auto iter2 = std::next(iter1); iter2 != mcParticles.end(); ++iter2)
            {
                const MCParticle *const pMC2{*iter2};
                if (pMC1->GetMomentum().GetMagnitude() < 1e-3 || pMC2->GetMomentum().GetMagnitude() < 1e-3)
                    continue;

                // Veto particles with too large an opening angle - suggest being conservative here, mainly want to avoid back-to-back cases
                const CartesianVector dir1{pMC1->GetMomentum().GetUnitVector()};
                const CartesianVector dir2{pMC2->GetMomentum().GetUnitVector()};
                const float costheta{dir1.GetDotProduct(dir2)};
                if (costheta < 0.866f)
                    continue;

                CaloHitList uHits2, vHits2, wHits2;
                this->CollectHitsByView(pMC2, uHits2, vHits2, wHits2);
                Eigen::MatrixXf uHitMatrix2; // Size will be set by FilterHits
                this->VectorizeAndFilterHits(uHits2, uVertex, m_vertexRadius, uHitMatrix2);
                Eigen::MatrixXf uHitMatrix1Filtered, uHitMatrix2Filtered;
                this->GetDifferenceAndFilterHits(uHitMatrix1, uHitMatrix2, m_distance, uHitMatrix1Filtered, uHitMatrix2Filtered);

                // Find shared hits by hashing one lot of filtered hits and storing in a set, then look for hashed hits from the other lot of filtered
                // hits in the set
                std::unordered_set<MatrixHit, MatrixHitHash> uHitSet;
                for (int i = 0; i < uHitMatrix1Filtered.rows(); ++i)
                {
                    uHitSet.insert({uHitMatrix1Filtered(i, 0), uHitMatrix1Filtered(i, 1)});
                }
                size_t sharedHits{0};
                for (int i = 0; i < uHitMatrix2Filtered.rows(); ++i)
                {
                    if (uHitSet.find({uHitMatrix2Filtered(i, 0), uHitMatrix2Filtered(i, 1)}) != uHitSet.end())
                        ++sharedHits;
                }
                const size_t uniqueHits1{uHitMatrix1Filtered.rows() - sharedHits};
                const size_t uniqueHits2{uHitMatrix2Filtered.rows() - sharedHits};
                const bool isCandidate{sharedHits > 2 || (uniqueHits1 >= 2 && uniqueHits2 >= 2)};
                if (isCandidate)
                    overlapCandidates[pMC1].insert(pMC2);
                if (m_visualise && isCandidate && uHitMatrix1Filtered.rows() > 0 && uHitMatrix2Filtered.rows() > 0)
                {
                    const CartesianVector pos(x, 0, u); (void)pos;
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits1, "1 (" + std::to_string(uHitMatrix1.rows()) + ")", BLACK));
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits2, "2 (" + std::to_string(uHitMatrix2.rows()) + ")", RED));
                    for (int i = 0; i < uHitMatrix1Filtered.rows(); ++i)
                    {
                        const CartesianVector pos1{uHitMatrix1Filtered(i, 0), 0, uHitMatrix1Filtered(i, 1)}; (void)pos1;
                        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos1, "1", BLACK, 1));
                    }
                    for (int i = 0; i < uHitMatrix2Filtered.rows(); ++i)
                    {
                        const CartesianVector pos2{uHitMatrix2Filtered(i, 0), 0, uHitMatrix2Filtered(i, 1)}; (void)pos2;
                        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos2, "2", RED, 1));
                    }

                    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "vtx", BLUE, 1));
                    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                }
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::AssessPfos(const MCToMCMap &overlapCandidates) const
{
    auto distanceSq = [](const CartesianVector &vertex, const CaloHit *const pCaloHit) -> float
    {
        const float dx{vertex.GetX() - pCaloHit->GetPositionVector().GetX()};
        const float dy{vertex.GetY() - pCaloHit->GetPositionVector().GetY()};
        const float dz{vertex.GetZ() - pCaloHit->GetPositionVector().GetZ()};
        return dx * dx + dy * dy + dz * dz;
    };

    // Look through the overlap candidates and identify the associated PFOs
    std::set<std::tuple<const MCParticle *, const Pfo *, const Pfo *>> pfoTuples;
    for (const auto &[pMC, pMCSet] : overlapCandidates)
    {
        if (pMCSet.empty())
            continue;
        const PfoSet &pfoSet{m_mcToPfoMap.at(pMC)};
        if (pfoSet.empty())
            continue;

        const CaloHitSet &caloHitSetMC1{m_mcToHitsMap.at(pMC)};
        CaloHitList caloHitListMC1;
        for (const CaloHit *const pCaloHit : caloHitSetMC1)
        {
            if (pCaloHit->GetHitType() == TPC_VIEW_U)
                caloHitListMC1.emplace_back(pCaloHit);
        }

        for (const MCParticle *const pMCOther : pMCSet)
        {
            const PfoSet &pfoSetOther{m_mcToPfoMap.at(pMCOther)};
            if (pfoSetOther.empty())
                continue;

            const CaloHitSet &caloHitSetMC2{m_mcToHitsMap.at(pMCOther)};
            CaloHitList caloHitListMC2;
            for (const CaloHit *const pCaloHit : caloHitSetMC2)
            {
                if (pCaloHit->GetHitType() == TPC_VIEW_U)
                    caloHitListMC2.emplace_back(pCaloHit);
            }

            for (const Pfo *const pPfo : pfoSet)
            {
                for (const Pfo *const pPfoOther : pfoSetOther)
                {
                    if (pPfo == pPfoOther)
                        continue;
                    pfoTuples.insert({pMC, pPfo, pPfoOther});
                    // Need to visualise these PFOs, or at least the hits sufficiently close to the common vertex
                    // to see what are inputs are. Seems like we're looking at quite a few PFOs given 2 MC

                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListMC1, "MC Hits 1", BLACK));
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListMC2, "MC Hits 2", RED));
                    CaloHitList caloHits1;
                    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHits1);
                    LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_U, caloHits1);
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits1, "PFO Hits 1", BLUE));
                    CaloHitList caloHits2;
                    LArPfoHelper::GetCaloHits(pPfoOther, TPC_VIEW_U, caloHits2);
                    LArPfoHelper::GetIsolatedCaloHits(pPfoOther, TPC_VIEW_U, caloHits2);
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits2, "PFO Hits 2", ORANGE));
                    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                }
            }
        }
    }

    const LArTransformationPlugin *pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (const auto &[pMC, pPfo1, pPfo2] : pfoTuples)
    {
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            if (view != TPC_VIEW_U)
                continue; // For now, only assess U view
            CaloHitList caloHits1;
            LArPfoHelper::GetCaloHits(pPfo1, view, caloHits1);
            LArPfoHelper::GetIsolatedCaloHits(pPfo1, view, caloHits1);
            if (caloHits1.empty())
                continue;
            CaloHitList caloHits2;
            LArPfoHelper::GetCaloHits(pPfo2, view, caloHits2);
            LArPfoHelper::GetIsolatedCaloHits(pPfo2, view, caloHits2);
            if (caloHits2.empty())
                continue;
            // Get the MC vertex and order the hits according to distance from the vertex
            float x{0}, c{0};
            LArVertexHelper::GetPositionProjection(pMC->GetVertex(), pTransform, view, x, c);
            const CartesianVector &vertex{x, 0, c};
            CaloHitVector caloHits1Filtered;
            std::copy_if(caloHits1.begin(), caloHits1.end(), std::back_inserter(caloHits1Filtered),
                [&](const CaloHit *hit) { return distanceSq(vertex, hit) < m_vertexRadius * m_vertexRadius; });
            if (caloHits1Filtered.empty())
                continue;
            CaloHitVector caloHits2Filtered;
            std::copy_if(caloHits2.begin(), caloHits2.end(), std::back_inserter(caloHits2Filtered),
                [&](const CaloHit *hit) { return distanceSq(vertex, hit) < m_vertexRadius * m_vertexRadius; });
            if (caloHits2Filtered.empty())
                continue;

            PcaResult pca1{this->PerformPca(caloHits1Filtered, vertex)};
            float meanDeviation1{std::accumulate(pca1.dT.begin(), pca1.dT.end(), 0.f) / pca1.dT.size()};
            if (meanDeviation1 > 1.f)
                continue;
            PcaResult pca2{this->PerformPca(caloHits2Filtered, vertex)};
            float meanDeviation2{std::accumulate(pca2.dT.begin(), pca2.dT.end(), 0.f) / pca2.dT.size()};
            if (meanDeviation2 > 1.f)
                continue;

            MahalanobisPairs mPairs;
            this->AlignPcaResults(pca1, pca2, mPairs);

            // Now we have a collection of hits to compare, assess the clusters to see if any hits should move between them
            std::vector<AssessmentResult> results;
            const StatusCode statusCode{this->AssessClusterAllocations(caloHits1Filtered, caloHits2Filtered, mPairs, results)};
            (void)statusCode; // Ignore the status code for now
            // Do something with the result
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::AssessClusterAllocations(const CaloHitVector &hits1, const CaloHitVector &hits2, const MahalanobisPairs &mPairs, std::vector<AssessmentResult> &results) const
{
    if (hits1.empty() || hits2.empty())
        return STATUS_CODE_SUCCESS;

    Eigen::MatrixXf hitMatrix1(hits1.size(), 2);
    LArEigenHelper::Vectorize(hits1, hitMatrix1);
    Eigen::MatrixXf hitMatrix2(hits2.size(), 2);
    LArEigenHelper::Vectorize(hits2, hitMatrix2);

    const LArTPC *const pTPC{this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second};
    const HitType view{hits1.front()->GetHitType()};
    const float pitch{view == TPC_VIEW_U ? pTPC->GetWirePitchU() : (view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW())};
    const float processVariance{m_processVarCoeff * pitch * pitch};
    const float measurementVariance{m_measurementVarCoeff * pitch * pitch};

    KalmanFilter2D kf1(KalmanFilter2D(m_delta, processVariance, measurementVariance, hitMatrix1.row(0).cast<double>().transpose()));
    //if (hits1.size() > 1)
    //    kf1.PredictAndUpdate(hitMatrix1.row(1).cast<double>().transpose());
    KalmanFilter2D kf2(KalmanFilter2D(m_delta, processVariance, measurementVariance, hitMatrix2.row(0).cast<double>().transpose()));
    //if (hits2.size() > 1)
    //    kf2.PredictAndUpdate(hitMatrix2.row(1).cast<double>().transpose());

    CaloHitList caloHitList1(hits1.begin(), hits1.end());
    CaloHitList caloHitList2(hits2.begin(), hits2.end());

    // Use the mPairs list to determine how to walk through the two sets of hits and whether or not to update the Kalman filters
    std::vector<AssessmentResult> results1;
    for (const auto &[i, i_reused, j, j_reused] : mPairs)
    {
        kf1.Predict();
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList1, "Hits 1", BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList2, "Hits 2", ORANGE));

        const KalmanFilter2D::PositionVector &hit1{hitMatrix1.row(i).cast<double>().transpose()};
        const KalmanFilter2D::PositionVector &hit2{hitMatrix2.row(j).cast<double>().transpose()};
        const CartesianVector pos1{static_cast<float>(hit1(0)), 0, static_cast<float>(hit1(1))};
        const CartesianVector pos2{static_cast<float>(hit2(0)), 0, static_cast<float>(hit2(1))};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos1, "1", BLACK, 1));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos2, "2", RED, 1));
        const double md11{kf1.GetMahalanobisDistance(hit1)};
        const double md12{kf1.GetMahalanobisDistance(hit2)};
        std::cout << "Tuple: (" << i << ", " << i_reused << ", " << j << ", " << j_reused << ")" << std::endl;
        std::cout << "MD11: " << md11 << ", MD12: " << md12 << std::endl;

        if (md12 < md11)
        {
            results1.emplace_back(AssessmentResult{hits2[j], 2, static_cast<int>(j), md12, 1});
            std::cout << "Red hit now part of blue cluster" << std::endl;
            kf1.Update(hit2);
        }
        if (!i_reused && md11 <= md12)
        {
            results1.emplace_back(AssessmentResult{hits1[i], 1, static_cast<int>(i), md11, 1});
            kf1.Update(hit1);
        }

        auto statePos1{kf1.GetTemporaryState().head<2>()};
        const CartesianVector pred1{static_cast<float>(statePos1(0)), 0, static_cast<float>(statePos1(1))};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pred1, "P", BLUE, 1));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    std::vector<AssessmentResult> results2;
    for (const auto &[i, i_reused, j, j_reused] : mPairs)
    {
        kf2.Predict();
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList1, "Hits 1", BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList2, "Hits 2", ORANGE));

        const KalmanFilter2D::PositionVector &hit1{hitMatrix1.row(i).cast<double>().transpose()};
        const KalmanFilter2D::PositionVector &hit2{hitMatrix2.row(j).cast<double>().transpose()};
        const CartesianVector pos1{static_cast<float>(hit1(0)), 0, static_cast<float>(hit1(1))};
        const CartesianVector pos2{static_cast<float>(hit2(0)), 0, static_cast<float>(hit2(1))};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos1, "1", BLACK, 1));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos2, "2", RED, 1));
        const double md21{kf2.GetMahalanobisDistance(hit1)};
        const double md22{kf2.GetMahalanobisDistance(hit2)};
        std::cout << "Tuple: (" << i << ", " << i_reused << ", " << j << ", " << j_reused << ")" << std::endl;
        std::cout << "MD21: " << md21 << ", MD22: " << md22 << std::endl;

        if (md21 < md22)
        {
            results2.emplace_back(AssessmentResult{hits1[i], 1, static_cast<int>(i), md21, 2});
            std::cout << "Black hit now part of orange cluster" << std::endl;
            kf2.Update(hit1);
        }
        if (!j_reused && md22 <= md21)
        {
            results2.emplace_back(AssessmentResult{hits2[j], 2, static_cast<int>(j), md22, 2});
            kf2.Update(hit2);
        }

        auto statePos1{kf2.GetTemporaryState().head<2>()};
        const CartesianVector pred1{static_cast<float>(statePos1(0)), 0, static_cast<float>(statePos1(1))};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pred1, "P", ORANGE, 1));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    // We've got a list of potential swaps now, but some of them may be claimed by both, pick the best
    for (const AssessmentResult &result1 : results1)
    {
        std::cout << "Assessing result 1: " << result1.pCaloHit << " Orig: " << result1.originalCluster << " MD: " << result1.md << " Suggested: " << result1.suggestedCluster << std::endl;
        bool found{false};
        for (const AssessmentResult &result2 : results2)
        {
            if (result1.pCaloHit == result2.pCaloHit)
            {
                if (result1.md < result2.md)
                {
                    results.emplace_back(result1);
                }
                else
                {
                    results.emplace_back(result2);
                }
                found = true;
                break;
            }
        }
        if (!found)
        {
            results.emplace_back(result1);
        }
    }
    for (const AssessmentResult &result2 : results2)
    {
        std::cout << "Assessing result 2: " << result2.pCaloHit << " Orig: " << result2.originalCluster << " MD: " << result2.md << " Suggested: " << result2.suggestedCluster << std::endl;
        bool found{false};
        for (const AssessmentResult &result1 : results1)
        {
            if (result2.pCaloHit == result1.pCaloHit)
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            results.emplace_back(result2);
        }
    }

    CaloHitList newList1, newList2;
    for (const AssessmentResult &result : results)
    {
        if (result.suggestedCluster == 1)
            newList1.emplace_back(result.pCaloHit);
        else if (result.suggestedCluster == 2)
            newList2.emplace_back(result.pCaloHit);
        std::cout << "Final result: " << result.pCaloHit << " Orig: " << result.originalCluster << " MD: " << result.md << " Suggested: " << result.suggestedCluster << std::endl;
    }

    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &newList1, "New 1", BLUE));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &newList2, "New 2", ORANGE));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapMonitoringAlgorithm::PcaResult TrackOverlapMonitoringAlgorithm::PerformPca(const CaloHitVector &hits, const CartesianVector &vertex) const
{
    PcaResult result;
    CartesianVector centroid{0.f, 0.f, 0.f};
    LArPcaHelper::EigenValues eigenValues{0.f, 0.f, 0.f};
    LArPcaHelper::EigenVectors eigenVectors;

    CartesianPointVector pcaPoints;
    result.centroid = vertex;
    for (const CaloHit *const pCaloHit : hits)
    {
        CartesianVector relativePos{pCaloHit->GetPositionVector() - vertex};
        pcaPoints.emplace_back(relativePos);
    }

    try
    {
        LArPcaHelper::RunPca(pcaPoints, centroid, eigenValues, eigenVectors);
    }
    catch (const StatusCodeException &)
    {
        std::cout << "PCA failed for " << hits.size() << " hits" << std::endl;
        return result;
    }

    result.principalAxis = eigenVectors.front();
    std::vector<std::pair<float, unsigned int>> projectedValues(hits.size());
    unsigned int negativeCount{0};
    for (unsigned int i = 0; i < hits.size(); ++i)
    {
        const CartesianVector &relativePos{pcaPoints[i]};
        projectedValues[i] = {relativePos.GetDotProduct(result.principalAxis), i};
        if (projectedValues[i].first < 0.f)
            ++negativeCount;
    }
    if (negativeCount > projectedValues.size() / 2)
    {
        // If more than half of the projected values are negative, flip the principal axis
        result.principalAxis = result.principalAxis * -1.f;
        for (auto &projectedValue : projectedValues)
        {
            projectedValue.first = -projectedValue.first;
        }
    }

    std::sort(projectedValues.begin(), projectedValues.end(),
        [](const std::pair<float, unsigned int> &lhs, const std::pair<float, unsigned int> &rhs) -> bool
        {
            return lhs.first < rhs.first;
        });
    std::cout << "Check sorted" << std::endl;
    for (const auto &[value, index] : projectedValues)
    {
        std::cout << "Hit " << index << ": " << value << std::endl;
    }

    for (const CartesianVector &relativePos : pcaPoints)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &relativePos, "PCA", GRAY, 1));
    }
    const CartesianVector a(0, 0, 0);
    const CartesianVector b(a + result.principalAxis * 5.f);
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &a, "vertex", BLACK, 2));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &a, &b, "PCA", ORANGE, 2, 1));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    result.sortedIndices.resize(projectedValues.size());
    result.dL.resize(projectedValues.size());
    result.dT.resize(projectedValues.size());
    for (unsigned int i = 0; i < projectedValues.size(); ++i)
    {
        result.sortedIndices[i] = projectedValues[i].second;
        result.dL[i] = projectedValues[i].first;
        // Compute the perpendicular distance using
        // || v - (v . u) * u || == sqrt(||v||^2 - (v . u)^2)
        const CartesianVector &v{pcaPoints[i]};
        const float vDotU{v.GetDotProduct(result.principalAxis)};
        const float vSq{v.GetMagnitudeSquared()};
        const float vDotUSq{vDotU * vDotU};
        if (vSq >= vDotUSq)
        {
            // Thise should generally be the case
            result.dT[i] = std::sqrt(vSq - vDotUSq);
        }
        else
        {
            // Can get here based on floating point precision issues
            const CartesianVector pointOnAxis{result.centroid + result.principalAxis * vDotU};
            result.dT[i] = (hits[i]->GetPositionVector() - pointOnAxis).GetMagnitude();
        }
    }
    result.succeeded = true;

    return result;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackOverlapMonitoringAlgorithm::AlignPcaResults(const PcaResult &pca1, const PcaResult &pca2, MahalanobisPairs &mPairs) const
{
    std::cout << "PCA 1: ";
    for (int index : pca1.sortedIndices)
    {
        std::cout << index << " ";
    }
    std::cout << std::endl;
    std::cout << "PCA 2: ";
    for (int index : pca2.sortedIndices)
    {
        std::cout << index << " ";
    }
    std::cout << std::endl;
    // For the Mahalanobis distance comparison to be meaningful, we can't walk through clusters at different rates
    // First, find the closest approach for each hit to a hit in the other cluster
    IntVector closestPartner1;
    for (size_t i = 0; i < pca1.dL.size(); ++i)
    {
        float minDistance{std::numeric_limits<float>::max()};
        int closestIndex{-1};
        for (size_t j = 0; j < pca2.dL.size(); ++j)
        {
            const float distance{std::abs(pca1.dL[i] - pca2.dL[j])};
            if (distance < minDistance)
            {
                minDistance = distance;
                closestIndex = static_cast<int>(j);
            }
        }
        closestPartner1.emplace_back(closestIndex);
    }
    IntVector closestPartner2;
    for (size_t i = 0; i < pca2.dL.size(); ++i)
    {
        float minDistance{std::numeric_limits<float>::max()};
        int closestIndex{-1};
        for (size_t j = 0; j < pca1.dL.size(); ++j)
        {
            const float distance{std::abs(pca2.dL[i] - pca1.dL[j])};
            if (distance < minDistance)
            {
                minDistance = distance;
                closestIndex = static_cast<int>(j);
            }
        }
        closestPartner2.emplace_back(closestIndex);
    }

    // Now create a list of candidate pairings
    std::vector<std::pair<int, int>> candidatePairs;
    for (size_t j = 0; j < closestPartner2.size(); ++j)
    {
        int i{closestPartner2[j]};
        candidatePairs.emplace_back(i, static_cast<int>(j));
    }
    for (size_t i = 0; i < closestPartner1.size(); ++i)
    {
        int j{closestPartner1[i]};
        candidatePairs.emplace_back(static_cast<int>(i), j);
    }

    // sort by i, then j
    std::sort(candidatePairs.begin(), candidatePairs.end());

    // Now construct an ordered set of pairs with duplciates removed, and track single index reusage
    std::unordered_set<std::pair<int, int>, PairHash> seenPairs;
    std::unordered_map<int, int> usage1, usage2;
    for (const auto &pair : candidatePairs)
    {
        int i{pair.first};
        int j{pair.second};
        if (seenPairs.insert(pair).second)
        {
            bool i_reused{usage1[i]++ > 0};
            bool j_reused{usage2[j]++ > 0};
            mPairs.emplace_back(std::make_tuple(i, i_reused, j, j_reused));
        }
    }

    for (const auto &[index1, isDuplicate1, index2, isDuplicate2] : mPairs)
    {
        std::cout << "Pair: " << index1 << " (" << isDuplicate1 << ") = " << pca1.dL[index1] << ", " << index2 << " (" << isDuplicate2 <<
            ") = " << pca2.dL[index2] << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackOverlapMonitoringAlgorithm::CollectHitsByView(const MCParticle *const pMC, CaloHitList &uHits, CaloHitList &vHits, CaloHitList &wHits) const
{
    const CaloHitSet &caloHits{m_mcToHitsMap.at(pMC)};
    for (const CaloHit *const pCaloHit : caloHits)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                uHits.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_V:
                vHits.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_W:
                wHits.emplace_back(pCaloHit);
                break;
            default:
                break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackOverlapMonitoringAlgorithm::VectorizeAndFilterHits(const pandora::CaloHitList &hits, const Eigen::RowVector2f &vertex, const float distance, Eigen::MatrixXf &filteredHits) const
{
    const float threshold{distance * distance};
    Eigen::MatrixXf allHitsMatrix(hits.size(), 2);
    LArEigenHelper::Vectorize(hits, allHitsMatrix);

    Eigen::MatrixXf diffsSq{(allHitsMatrix.rowwise() - vertex).rowwise().squaredNorm()};
    std::vector<int> keepIndices;
    for (int i = 0; i < diffsSq.rows(); ++i)
    {
        if (diffsSq(i, 0) < threshold)
            keepIndices.emplace_back(i);
    }
    filteredHits.resize(keepIndices.size(), allHitsMatrix.cols());
    for (int i = 0; i < static_cast<int>(keepIndices.size()); ++i)
        filteredHits.row(i) = allHitsMatrix.row(keepIndices[i]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackOverlapMonitoringAlgorithm::GetDifferenceAndFilterHits(const Eigen::MatrixXf &hits1, const Eigen::MatrixXf &hits2, const float distance, Eigen::MatrixXf &filteredHits1, Eigen::MatrixXf &filteredHits2) const
{
    const float threshold{distance * distance};
    Eigen::VectorXf sqNorm1(hits1.rowwise().squaredNorm());
    Eigen::VectorXf sqNorm2(hits2.rowwise().squaredNorm());
    Eigen::MatrixXf distances = sqNorm1.replicate(1, hits2.rows()) + sqNorm2.transpose().replicate(hits1.rows(), 1) - 2.f * (hits1 * hits2.transpose());

    std::unordered_set<int> keepIndices1, keepIndices2;
    for (int i = 0; i < distances.rows(); ++i)
    {
        for (int j = 0; j < distances.cols(); ++j)
        {
            if (distances(i, j) < threshold)
            {
                keepIndices1.insert(i);
                keepIndices2.insert(j);
            }
        }
    }

    filteredHits1.resize(keepIndices1.size(), hits1.cols());
    auto iter1 = keepIndices1.begin();
    for (int i = 0; i < static_cast<int>(keepIndices1.size()); ++i, ++iter1)
        filteredHits1.row(i) = hits1.row(*iter1);
    filteredHits2.resize(keepIndices2.size(), hits2.cols());
    auto iter2 = keepIndices2.begin();
    for (int i = 0; i < static_cast<int>(keepIndices2.size()); ++i, ++iter2)
        filteredHits2.row(i) = hits2.row(*iter2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));

    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRadius", m_vertexRadius));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OverlapDistance", m_distance));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Delta", m_delta));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProcessVarianceCoefficient", m_processVarCoeff));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MeasurementVarianceCoefficient", m_measurementVarCoeff));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
