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

    CaloHitList refHitsU, refHitsV, refHitsW, refHits;
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
                    refHits.emplace_back(pCaloHit);
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
                    refHits.emplace_back(pCaloHit);
                    refPositions3D.emplace_back(position3D);
                }
            }
        }
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

        for (const CaloHit *const pCaloHit : refHits)
        {
            const LArPlaneContextObject::HitTriplet *const pHitTriplet(m_planeContextObject->GetHitTriplet(pCaloHit));
            if (pHitTriplet && pHitTriplet->m_chi2 < 0.1f)
            {
                const CartesianVector &position(pHitTriplet->m_wHit->GetPositionVector());
                try
                {
                    const MCParticle *const pMC(MCParticleHelper::GetMainMCParticle(pHitTriplet->m_wHit));
                    if (pMC->GetParticleId() == 22)
                        continue;
                    //std::cout << "Found triplet: " << pMC->GetParticleId() << " (" << position.GetX() << ", " << position.GetY() << ", "
                    //    << position.GetZ() << ") " << "chi2: " << pHitTriplet->m_chi2 << std::endl;
                }
                catch (StatusCodeException &)
                {
                }
            }
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

StatusCode CaloMatchTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PlaneSolverContextName",
        m_planeSolverContextName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
