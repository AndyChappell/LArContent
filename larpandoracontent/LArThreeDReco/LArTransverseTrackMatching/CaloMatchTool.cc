/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CaloMatchTool.cc
 *
 *  @brief  Implementation of the calorimetric matching tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArDiscreteProbabilityHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CaloMatchTool.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

using namespace pandora;

namespace lar_content
{

CaloMatchTool::CaloMatchTool() :
    m_planeContextObject(nullptr),
    m_planeSolverContextName("PlaneContext")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CaloMatchTool::Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    m_planeContextObject = dynamic_cast<const LArPlaneContextObject *>(PandoraContentApi::GetEventContextObject(*pAlgorithm, m_planeSolverContextName));

    if (!m_planeContextObject)
        return false;

    ClusterVector keyClusters;
    overlapTensor.GetSortedKeyClusters(keyClusters);

    for (const Cluster *const pKeyCluster : keyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        for (const auto &element : elementList)
        {
            const Cluster *const pClusterU(element.GetClusterU());
            const Cluster *const pClusterV(element.GetClusterV());
            const Cluster *const pClusterW(element.GetClusterW());

            Trajectory trajectory;
            if (!this->Construct3DTrajectoryTemplates(pClusterU, pClusterV, pClusterW, trajectory))
                continue;

            CaloMatchResult caloResult;
            if (!this->CompareChargeProfiles(trajectory, caloResult))
                continue;

            std::cout << "Calo comparison: " << "nBins=" << caloResult.m_nBins << " rUV=" << caloResult.m_uv.m_correlation << " pUV="
                << caloResult.m_uv.m_pValue << " rUW=" << caloResult.m_uw.m_correlation << " pUW=" << caloResult.m_uw.m_pValue << " rVW="
                << caloResult.m_vw.m_correlation << " pVW=" << caloResult.m_vw.m_pValue << " goodMatch=" << caloResult.m_isGoodMatch << std::endl;
        }
    }

    return false; // Placeholder
    //ProtoParticleVector protoParticleVector;
    //this->FindLongTracks(overlapTensor, protoParticleVector);

    //const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    //return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CaloMatchTool::Construct3DTrajectoryTemplates(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
    const pandora::Cluster *const pClusterW, Trajectory &trajectory) const
{
    CaloHitList refHitsU, refHitsV, refHitsW;
    CartesianPointVector refPositions3D;
    LArClusterHelper::GetAllHits(pClusterU, refHitsU);
    LArClusterHelper::GetAllHits(pClusterV, refHitsV);
    LArClusterHelper::GetAllHits(pClusterW, refHitsW);
    if (refHitsU.empty() || refHitsV.empty() || refHitsW.empty())
        return false;
    
    this->Estimate3DTrajectory(refHitsU, refHitsV, refHitsW, refPositions3D);
    CartesianPointVector trajectoryTemplate;
    const unsigned int nKeyPoints = 10;
    this->SparsifyTrajectory(refPositions3D, nKeyPoints, trajectoryTemplate, refHitsW);

    if (trajectoryTemplate.size() < 2)
        return false;

    const float binSize = 2.f;
    this->ProjectTrajectory(trajectoryTemplate, refHitsW, binSize, trajectory);
    this->ProjectTrajectory(trajectoryTemplate, refHitsU, binSize, trajectory);
    this->ProjectTrajectory(trajectoryTemplate, refHitsV, binSize, trajectory);

    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &refHitsU, "U Hits", RED));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &refHitsV, "V Hits", GREEN));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &refHitsW, "W Hits", BLUE));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return !trajectory.empty();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloMatchTool::Estimate3DTrajectory(const CaloHitList &refHitsU, const CaloHitList &refHitsV, const CaloHitList &refHitsW,
    CartesianPointVector &refPositions3D) const
{
    for (const CaloHit *const pCaloHit : refHitsW)
    {
        const LArPlaneContextObject::HitTriplet *const pHitTriplet(m_planeContextObject->GetHitTriplet(pCaloHit));
        if (pHitTriplet && pHitTriplet->m_chi2 < 1.f)
        {
            const CaloHit *pHitU(nullptr), *pHitV(nullptr);
            if (std::find(refHitsU.begin(), refHitsU.end(), pHitTriplet->m_uHit) != refHitsU.end())
                pHitU = pHitTriplet->m_uHit;
            if (std::find(refHitsV.begin(), refHitsV.end(), pHitTriplet->m_vHit) != refHitsV.end())
                pHitV = pHitTriplet->m_vHit;
            if (pHitTriplet->m_nHits == 3 && pHitU && pHitV)
            {
                CartesianVector position3D(0.f, 0.f, 0.f);
                float chi2(0.f);
                LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W,
                    pHitU->GetPositionVector(), pHitV->GetPositionVector(), pCaloHit->GetPositionVector(), position3D, chi2);
                refPositions3D.emplace_back(position3D);
            }
            else
            {
                const CaloHit *pOtherHit(pHitU ? pHitU : pHitV);
                if (!pOtherHit)
                    continue;
                CartesianVector position3D(0.f, 0.f, 0.f);
                float chi2(0.f);
                LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_W, pOtherHit->GetHitType(),
                    pHitTriplet->m_wHit->GetPositionVector(), pOtherHit->GetPositionVector(), position3D, chi2);
                refPositions3D.emplace_back(position3D);
            }
        }
    }
}

void CaloMatchTool::SparsifyTrajectory(const CartesianPointVector &refPositions3D, const unsigned int nKeyPoints,
    CartesianPointVector &trajectoryTemplate, const CaloHitList &refHitsW) const
{
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(refPositions3D, centroid, eigenValues, eigenVectors);
    LArPcaHelper::ProjectedPointVector projectedPoints(LArPcaHelper::GetOrderedPoints(refPositions3D));

    // Order using the 3D principal axis projected into the W view,
    const LArTransformationPlugin *const pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};

    const auto getW = [&](const CartesianVector &p) -> float
    {
        return static_cast<float>(pTransform->YZtoW(p.GetY(), p.GetZ()));
    };

    const float centroidW = getW(centroid);
    const float axisX = eigenVectors[0].GetX();
    const float axisW = getW(centroid + eigenVectors[0]) - centroidW;

    for (auto &[l, position3D] : projectedPoints)
        l = (position3D.GetX() - centroid.GetX()) * axisX + (getW(position3D) - centroidW) * axisW;

    std::stable_sort(projectedPoints.begin(), projectedPoints.end(),
        [](const auto &a, const auto &b)
        {
            return a.first < b.first;
        });

    for (const auto &[l, position3D] : projectedPoints)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position3D, "Pos: " + std::to_string(l), BLUE, 2));
    }
    CartesianVector start(centroid - eigenVectors[0] * 20);
    CartesianVector end(centroid + eigenVectors[0] * 20);
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &refHitsW, "2D Hits", RED));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &end, "Fit", BLACK, 3, 1));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    // Extract a spread of points along the ordered trajectory
    // ATTN: Need to think about forward/backward directionality depending upon the flipping of clusters (for now, it's arbitrary)
    if (projectedPoints.size() > nKeyPoints && nKeyPoints >= 2)
    {
        const float tMin = projectedPoints.front().first;
        const float tMax = projectedPoints.back().first;

        size_t firstAllowed(0);

        for (unsigned int i = 0; i < nKeyPoints; ++i)
        {
            const float target = tMin + (tMax - tMin) * static_cast<float>(i) / static_cast<float>(nKeyPoints - 1);

            // Leave sufficient points for all remaining selections
            const size_t lastAllowed = projectedPoints.size() - (nKeyPoints - i);

            size_t bestIndex(firstAllowed);
            float bestDelta(std::numeric_limits<float>::max());

            for (size_t j = firstAllowed; j <= lastAllowed; ++j)
            {
                const float delta = std::abs(projectedPoints[j].first - target);
                if (delta < bestDelta)
                {
                    bestDelta = delta;
                    bestIndex = j;
                }
            }

            trajectoryTemplate.emplace_back(projectedPoints[bestIndex].second);

            // Future selections must use subsequent points
            firstAllowed = bestIndex + 1;
        }
    }
    else
    {
        for (const auto &[l, position3D] : projectedPoints)
            trajectoryTemplate.emplace_back(position3D);
    }

    for (const CartesianVector &position3D : trajectoryTemplate)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position3D, "Template", BLUE, 2));
    }
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &refHitsW, "2D Hits", RED));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloMatchTool::ProjectTrajectory(const CartesianPointVector &trajectoryTemplate, const CaloHitList &refHits, const float binSize, Trajectory &trajectory) const
{
    if (refHits.empty() || trajectoryTemplate.size() < 2 || binSize <= 0.f)
        return;

    const HitType view = refHits.front()->GetHitType();
    if (view != TPC_VIEW_U && view != TPC_VIEW_V && view != TPC_VIEW_W)
        return;

    const LArTransformationPlugin *const pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};

    // Project a 3D point into the requested 2D view.
    const auto project = [&](const CartesianVector &p) -> CartesianVector
    {
        float coordinate(0.f);
        switch (view)
        {
            case TPC_VIEW_U:
                coordinate = static_cast<float>(pTransform->YZtoU(p.GetY(), p.GetZ()));
                break;
            case TPC_VIEW_V:
                coordinate = static_cast<float>(pTransform->YZtoV(p.GetY(), p.GetZ()));
                break;
            case TPC_VIEW_W:
                coordinate = static_cast<float>(pTransform->YZtoW(p.GetY(), p.GetZ()));
                break;
            default:
                break;
        }

        return CartesianVector(p.GetX(), 0.f, coordinate);
    };

    // Piece3D describes one straight portion of the trajectory. Multiple pieces can belong to the same bin (e.g. around initial template point).
    // BinGeometry describes the complete bin, whose endpoints may span one or more changes in trajectory direction.
    struct Piece3D
    {
        CartesianVector m_start;
        CartesianVector m_end;
        int m_binIndex;
    };

    struct BinGeometry
    {
        CartesianVector m_start;
        CartesianVector m_end;
        double m_pathLength;
    };

    std::vector<Piece3D> pieces3D;
    std::vector<BinGeometry> binDefinitions;

    const double targetLength = static_cast<double>(binSize);

    const double tolerance = 1.e-3;
    double distanceInBin(0.);
    int binIndex(0);

    // Build the trajectory pieces and bin geometry in 3D (view-independent).
    for (size_t i = 1; i < trajectoryTemplate.size(); ++i)
    {
        const CartesianVector &start3D = trajectoryTemplate[i - 1];
        const CartesianVector &end3D = trajectoryTemplate[i];
        const double dx = static_cast<double>(end3D.GetX()) - start3D.GetX();
        const double dy = static_cast<double>(end3D.GetY()) - start3D.GetY();
        const double dz = static_cast<double>(end3D.GetZ()) - start3D.GetZ();
        const double segmentLength = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (segmentLength <= 0.)
            continue;

        const auto interpolate = [&](const double t) -> CartesianVector
        {
            return CartesianVector(static_cast<float>(start3D.GetX() + dx * t), static_cast<float>(start3D.GetY() + dy * t),
                static_cast<float>(start3D.GetZ() + dz * t));
        };

        double distanceInSegment(0.);
        while (distanceInSegment < segmentLength)
        {
            const double segmentRemaining = segmentLength - distanceInSegment;
            const double binRemaining = targetLength - distanceInBin;

            // If the template vertex is very close to the next bin boundary, snap the boundary onto that vertex.
            const bool atBoundary = std::abs(segmentRemaining - binRemaining) <= tolerance;
            const double step = atBoundary ? segmentRemaining : std::min(segmentRemaining, binRemaining);
            const bool completesBin = (segmentRemaining >= binRemaining - tolerance);

            // Interpolate the piece endpoints along the current straight section of the original 3D trajectory.
            const double t0 = distanceInSegment / segmentLength;
            const double t1 = (distanceInSegment + step) / segmentLength;

            const CartesianVector p0(interpolate(t0));
            const CartesianVector p1(interpolate(t1));
            pieces3D.push_back({p0, p1, binIndex});

            // Construct or extend the corresponding complete 3D bin.
            if (static_cast<size_t>(binIndex) == binDefinitions.size())
            {
                binDefinitions.push_back({p0, p1, step});
            }
            else
            {
                BinGeometry &bin = binDefinitions.at(binIndex);
                bin.m_end = p1;
                bin.m_pathLength += step;
            }

            if (step >= segmentRemaining)
                distanceInSegment = segmentLength;
            else
                distanceInSegment += step;

            if (completesBin)
            {
                // The next piece belongs to the next bin.
                ++binIndex;
                distanceInBin = 0.;
            }
            else
            {
                // Continue this bin across the template vertex.
                distanceInBin += step;
            }
        }
    }

    // Initialise the shared trajectory on the first call. On subsequent calls, ensure the supplied 3D template and bin size generate matching bin geometry.
    if (trajectory.empty())
    {
        for (size_t i = 0; i < binDefinitions.size(); ++i)
        {
            const BinGeometry &bin = binDefinitions[i];
            trajectory.emplace(static_cast<int>(i), std::make_unique<TrajectorySegment>(bin.m_start, bin.m_end, static_cast<float>(bin.m_pathLength)));
        }
    }
    else
    {
        if (trajectory.size() != binDefinitions.size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        constexpr float geometryTolerance = 1.e-3f;

        for (size_t i = 0; i < binDefinitions.size(); ++i)
        {
            const auto iter = trajectory.find(static_cast<int>(i));

            if (iter == trajectory.end())
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            const TrajectorySegment &existing = *iter->second;
            const BinGeometry &expected = binDefinitions[i];

            const float startDifference = (existing.m_start3D - expected.m_start).GetMagnitude();
            const float endDifference = (existing.m_end3D - expected.m_end).GetMagnitude();

            const double lengthDifference = std::abs(static_cast<double>(existing.m_pathLength3D) - expected.m_pathLength);

            if (startDifference > geometryTolerance || endDifference > geometryTolerance || lengthDifference > geometryTolerance)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }

    // Project the individual 3D pieces into the requested view.
    struct Piece2D
    {
        CartesianVector m_start;
        CartesianVector m_end;
        int m_binIndex;
    };

    std::vector<Piece2D> pieces2D;
    pieces2D.reserve(pieces3D.size());

    for (const Piece3D &piece : pieces3D)
        pieces2D.push_back({project(piece.m_start), project(piece.m_end), piece.m_binIndex});

    // Reset only the charge belonging to this view.
    for (auto &[index, pSegment] : trajectory)
    {
        switch (view)
        {
            case TPC_VIEW_U:
                pSegment->m_chargeU = 0.f;
                break;
            case TPC_VIEW_V:
                pSegment->m_chargeV = 0.f;
                break;
            case TPC_VIEW_W:
                pSegment->m_chargeW = 0.f;
                break;
            default:
                break;
        }
    }

    // Associate each hit with the nearest projected piece.
    for (const CaloHit *const pCaloHit : refHits)
    {
        if (pCaloHit->GetHitType() != view)
            continue;

        const CartesianVector &position(pCaloHit->GetPositionVector());
        float bestDistanceSquared(std::numeric_limits<float>::max());
        int bestBinIndex(-1);
        for (const Piece2D &piece : pieces2D)
        {
            const CartesianVector delta(piece.m_end - piece.m_start);
            const float lengthSquared = delta.GetMagnitudeSquared();
            float t(0.f);

            if (lengthSquared > 1.e-10f)
            {
                // Project onto the finite straight piece.
                t = (position - piece.m_start).GetDotProduct(delta) / lengthSquared;
                t = std::clamp(t, 0.f, 1.f);
            }

            const CartesianVector closestPoint(piece.m_start + delta * t);
            const float distanceSquared = (position - closestPoint).GetMagnitudeSquared();

            if (distanceSquared < bestDistanceSquared)
            {
                bestDistanceSquared = distanceSquared;
                bestBinIndex = piece.m_binIndex;
            }
        }

        if (bestBinIndex < 0)
            continue;

        // The selected bin is shared across all three views.
        TrajectorySegment &bin = *trajectory.at(bestBinIndex);
        const float charge = pCaloHit->GetInputEnergy();
        switch (view)
        {
            case TPC_VIEW_U:
                bin.m_chargeU += charge;
                break;
            case TPC_VIEW_V:
                bin.m_chargeV += charge;
                break;
            case TPC_VIEW_W:
                bin.m_chargeW += charge;
                break;
            default:
                break;
        }
    }

    // Visualise the projected boundaries of the complete 3D bins.
    bool firstSegment(true);

    for (const auto &[index, pSegment] : trajectory)
    {
        if (firstSegment)
        {
            const CartesianVector position2D(project(pSegment->m_start3D));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position2D, "traj", MAGENTA, 1));
            firstSegment = false;
        }

        const CartesianVector position2D(project(pSegment->m_end3D));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position2D, "traj", MAGENTA, 1));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------


bool CaloMatchTool::CompareChargeProfiles(const Trajectory &trajectory, CaloMatchResult &result) const
{
    const unsigned int minBins = 2;
    const unsigned int nPermutations = 1000;
    const float minCorrelation = 0.5f;
    const float maxPValue = 0.05f;

    DiscreteProbabilityVector::AllFloatInputData inputU, inputV, inputW;

    inputU.reserve(trajectory.size());
    inputV.reserve(trajectory.size());
    inputW.reserve(trajectory.size());

    float pathLength = 0.f;
    float totalU = 0.f, totalV = 0.f, totalW = 0.f;

    for (const auto &[index, pSegment] : trajectory)
    {
        const float length = pSegment->m_pathLength3D;
        if (length <= 0.f)
            return false;

        // Only include bins where at least one view has charge
        if (pSegment->m_chargeU + pSegment->m_chargeV + pSegment->m_chargeW > 0.f)
        {
            inputU.emplace_back(pathLength, pSegment->m_chargeU);
            inputV.emplace_back(pathLength, pSegment->m_chargeV);
            inputW.emplace_back(pathLength, pSegment->m_chargeW);

            totalU += pSegment->m_chargeU;
            totalV += pSegment->m_chargeV;
            totalW += pSegment->m_chargeW;
        }

        pathLength += length;
    }

    result.m_nBins = static_cast<unsigned int>(inputU.size());
    if (result.m_nBins < minBins)
        return false;

    // The probability-vector constructor requires a positive normalisation for every view.
    const float epsilon = std::numeric_limits<float>::epsilon();
    if (totalU <= epsilon || totalV <= epsilon || totalW <= epsilon)
        return false;

    try
    {
        const DiscreteProbabilityVector profileU(inputU, pathLength, false);
        const DiscreteProbabilityVector profileV(inputV, pathLength, false);
        const DiscreteProbabilityVector profileW(inputW, pathLength, false);

        std::mt19937 randomNumberGenerator(42);
        const auto comparePair = [&](const DiscreteProbabilityVector &a, const DiscreteProbabilityVector &b) -> CaloPairResult
        {
            CaloPairResult pair;
            pair.m_correlation = LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(a, b);
            pair.m_pValue = LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromPermutationTest(a, b, randomNumberGenerator, nPermutations);
            return pair;
        };

        result.m_uv = comparePair(profileU, profileV);
        result.m_uw = comparePair(profileU, profileW);
        result.m_vw = comparePair(profileV, profileW);

        const float minR = std::min({result.m_uv.m_correlation, result.m_uw.m_correlation, result.m_vw.m_correlation});
        const float maxP = std::max({result.m_uv.m_pValue, result.m_uw.m_pValue, result.m_vw.m_pValue});
        result.m_isGoodMatch = std::isfinite(minR) && std::isfinite(maxP) && minR >= minCorrelation && maxP <= maxPValue;

        return true;
    }
    catch (const pandora::StatusCodeException &)
    {
        return false;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloMatchTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PlaneSolverContextName",
        m_planeSolverContextName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
