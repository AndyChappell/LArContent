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
    const LArTransformationPlugin *const pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};

    CaloHitList refHitsU, refHitsV, refHitsW;
    CartesianPointVector refPositions3D;
    LArClusterHelper::GetAllHits(pClusterU, refHitsU);
    LArClusterHelper::GetAllHits(pClusterV, refHitsV);
    LArClusterHelper::GetAllHits(pClusterW, refHitsW);
    std::cout << "Ref hits: " << refHitsU.size() << ", " << refHitsV.size() << ", " << refHitsW.size() << std::endl;
    if (!refHitsW.empty())
    {
        // Find all the hit triplets that are fully consistent with the three clusters, and have a good chi2.
        // Then, apply a PCA fit to those points to get a baseline direction
        // Then, filter to only those hits with excellent chi2 and try to select a roughly even spread of the hits along the baseline
        // trajectory. These will form the trajectory template
        this->Estimate3DTrajectory(refHitsU, refHitsV, refHitsW, refPositions3D);
        // PCA fit the refPositions3D here
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
        CartesianPointVector trajectoryTemplate;
        const unsigned int nKeyPoints = 10;
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

        // Project the template into each 2D view
        CartesianPointVector trajectoryTemplateW;
        for (const CartesianVector &position3D : trajectoryTemplate)
        {
            const CartesianVector projection(position3D.GetX(), 0.f, static_cast<float>(pTransform->YZtoW(position3D.GetY(), position3D.GetZ())));
            trajectoryTemplateW.emplace_back(projection);
        }

        // Build the trajectory segment objects and populate with the relevant hits in each view.
        Trajectory trajectory;
        for (const CaloHit *const pCaloHit : refHitsW)
        {
            const CartesianVector &position(pCaloHit->GetPositionVector());
            float bestDistanceSquared(std::numeric_limits<float>::max());
            float bestLongitudinal(std::numeric_limits<float>::max());
            int bestSegmentIndex(-1);
            for (size_t i = 0; i < trajectoryTemplateW.size() - 1; ++i)
            {
                const bool firstSegment(i == 0), lastSegment(i == trajectoryTemplateW.size() - 2);
                const CartesianVector &start(trajectoryTemplateW[i]);
                const CartesianVector &end(trajectoryTemplateW[i + 1]);
                const CartesianVector axis = (end - start).GetUnitVector();
                const float limitLow = 0.f;
                const float limitHigh = (end - start).GetMagnitude();
                const float projection = (position - start).GetDotProduct(axis);
                const CartesianVector closestPoint = start + axis * projection;
                if ((projection >= limitLow || firstSegment) && (projection < limitHigh || lastSegment))
                {
                    const float distanceSquared = (position - closestPoint).GetMagnitudeSquared();
                    if (distanceSquared < bestDistanceSquared)
                    {
                        bestDistanceSquared = distanceSquared;
                        bestLongitudinal = projection;
                        bestSegmentIndex = i;
                    }
                }
            }

            if (bestSegmentIndex >= 0)
            {
                if (trajectory.find(bestSegmentIndex) == trajectory.end())
                {
                    trajectory[bestSegmentIndex] = std::make_unique<TrajectorySegment>(trajectoryTemplateW[bestSegmentIndex],
                        trajectoryTemplateW[bestSegmentIndex + 1], ChargeProfile(), ChargeProfile(), ChargeProfile());
                }
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                        trajectory[bestSegmentIndex]->m_chargeProfileU.emplace_back(bestLongitudinal, pCaloHit->GetInputEnergy());
                        break;
                    case TPC_VIEW_V:
                        trajectory[bestSegmentIndex]->m_chargeProfileV.emplace_back(bestLongitudinal, pCaloHit->GetInputEnergy());
                        break;
                    case TPC_VIEW_W:
                        trajectory[bestSegmentIndex]->m_chargeProfileW.emplace_back(bestLongitudinal, pCaloHit->GetInputEnergy());
                        break;
                    default:
                        break;
                }
            }
        }
        for (const auto &[index, pSegment] : trajectory)
        {
            if (!pSegment)
                continue;
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &refHitsW, "2D Hits", RED));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &pSegment->m_start, &pSegment->m_end, "Template", BLACK, 3, 1));
            for (const auto &[l, adc] : pSegment->m_chargeProfileW)
            {
                CartesianVector pos(pSegment->m_start + (pSegment->m_end - pSegment->m_start).GetUnitVector() * l);
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, std::to_string(index) + " L: " + std::to_string(l)
                    + " Q: " + std::to_string(adc), GREEN, 2));
            }
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
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

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloMatchTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PlaneSolverContextName",
        m_planeSolverContextName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
