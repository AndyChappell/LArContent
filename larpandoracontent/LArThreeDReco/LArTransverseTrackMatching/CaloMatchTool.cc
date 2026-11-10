/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CaloMatchTool.cc
 *
 *  @brief  Implementation of the calorimetric matching tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CaloMatchTool.h"

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

            pandora::CartesianPointVector trajectory;
            if (this->Construct3DTrajectoryTemplates(pClusterU, pClusterV, pClusterW, trajectory))
            {
                //std::cout << "Found trajectory" << std::endl;
            }
            else
            {
                //std::cout << "Failed to find trajectory" << std::endl;
            }
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
    const pandora::Cluster *const pClusterW, pandora::CartesianPointVector &trajectory) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    CaloHitList refHitsU, refHitsV, refHitsW;
    CartesianPointVector refPositions3D;
    LArClusterHelper::GetAllHits(pClusterU, refHitsU);
    LArClusterHelper::GetAllHits(pClusterV, refHitsV);
    LArClusterHelper::GetAllHits(pClusterW, refHitsW);
    std::cout << "Ref hits: " << refHitsU.size() << ", " << refHitsV.size() << ", " << refHitsW.size() << std::endl;
    if (!refHitsW.empty())
    {
        this->Estimate3DTrajectory(refHitsU, refHitsV, refHitsW, refPositions3D);
        CartesianPointVector trajectoryTemplate;
        const unsigned int nKeyPoints = 10;
        this->SparsifyTrajectory(refPositions3D, nKeyPoints, trajectoryTemplate, refHitsW);   // refHitsW is temporary for debugging visualisations

        Trajectory trajectory;
        this->ProjectTrajectory(trajectoryTemplate, refHitsW, 2.f, trajectory);
        // Need to add U and V as well
    }
    // Get the hits from each cluster.
    // Find the full triplet matches with an excellent chi2.
    // Figure out the ordering of those triplets.
    // Pick a sensible spreads of those triplets to make a trajectory template.
    // We want the minimum number of points that define reasonably straight segments.
    // Then aim to put together geometric and calorimetric consistency between the clusters based on the trajectory template.

    return false;
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
    if (projectedPoints.size() > nKeyPoints)
    {
        const float tMin = projectedPoints.front().first;
        const float tMax = projectedPoints.back().first;
        auto it = projectedPoints.begin();
        for (unsigned int i = 0; i < nKeyPoints; ++i)
        {
            float bestDelta(std::numeric_limits<float>::max());
            const float target = tMin + (tMax - tMin) * static_cast<float>(i) / static_cast<float>(nKeyPoints - 1);
            while (true)
            {
                const float delta = std::abs(it->first - target);
                if (delta < bestDelta)
                {
                    bestDelta = delta;
                    ++it;
                }
                else
                {
                    --it;
                    if (std::find(trajectoryTemplate.begin(), trajectoryTemplate.end(), it->second) == trajectoryTemplate.end())
                        trajectoryTemplate.emplace_back(it->second);
                    break;
                }
                if (bestDelta == std::numeric_limits<float>::max() && it == projectedPoints.end())
                    break;
            }
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

StatusCode CaloMatchTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PlaneSolverContextName",
        m_planeSolverContextName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
