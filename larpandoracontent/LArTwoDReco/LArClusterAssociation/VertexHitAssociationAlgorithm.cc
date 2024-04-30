/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/VertexHitAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/VertexHitAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexHitAssociationAlgorithm::VertexHitAssociationAlgorithm() :
    m_minClusterLayers(4),
    m_maxGapLayers(7),
    m_fitLayers(30),
    m_maxGapDistanceSquared(10.f),
    m_minCosRelativeAngle(0.985f),
    m_maxTransverseDisplacement(2.f),
    m_maxLongitudinalDisplacement(2.f),
    m_hitSizeZ(0.3f),
    m_hitSizeX(0.5f),
    m_view(TPC_3D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexHitAssociationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    ClusterList selectedClusters;
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetOrderedCaloHitList().size() > 5)
            selectedClusters.emplace_back(pCluster);
    }

    // Get the hits and the clusters
    // Project the cluster end direction to see how close it gets to various hits
    // Keep plausible candiate matches
    // Consider groups of hits (in close proximity and "matched" to some fit trajectory
    // Decide which association is best

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    std::vector<LArPointingCluster> pointingClusterVector;
    std::map<const Cluster *, CaloHitList> innerClusterToCandidateMap, outerClusterToCandidateMap;
    std::map<const Cluster *, int> clusterToIndexMap;
    std::map<const CaloHit *, int> hitAssocMap;
    for (const Cluster *const pCluster : selectedClusters)
    {
        HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        ClusterList clusterList;
        clusterList.emplace_back(pCluster);
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "cluster", BLACK, false));
        const LArPointingCluster pointingCluster(pCluster, 10, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
        const CartesianVector &a{pointingCluster.GetInnerVertex().GetPosition()};
        const CartesianVector &b{a - pointingCluster.GetInnerVertex().GetDirection() * 10};
        const CartesianVector &c{pointingCluster.GetOuterVertex().GetPosition()};
        const CartesianVector &d{c - pointingCluster.GetOuterVertex().GetDirection() * 10};
        const CartesianVector &innerDir{pointingCluster.GetInnerVertex().GetDirection() * -1.f};
        const CartesianVector &outerDir{pointingCluster.GetOuterVertex().GetDirection() * -1.f};
        pointingClusterVector.emplace_back(pointingCluster);
        clusterToIndexMap[pCluster] = pointingClusterVector.size() - 1;

        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &a, &b, "ab", BLUE, 1, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &c, &d, "cd", RED, 1, 1));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            const CartesianVector &vecInner{pos - pointingCluster.GetInnerVertex().GetPosition()};
            if (vecInner.GetMagnitude() <= 10.f)
            {
                const float innerDot{innerDir.GetDotProduct(vecInner)};
                if (innerDot >= -0.5f && innerDot <= 10.f)
                {
                    const float innerCross{innerDir.GetCrossProduct(vecInner).GetMagnitude()};
                    if (innerCross < 1.f)
                    {
                        innerClusterToCandidateMap[pCluster].emplace_back(pCaloHit);
                        if (hitAssocMap.find(pCaloHit) == hitAssocMap.end())
                            hitAssocMap[pCaloHit] = 1;
                        else
                            ++hitAssocMap[pCaloHit];
                    }
                }
                //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pCaloHit->GetPositionVector(), "hit", BLACK, 1));
                //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }

            const CartesianVector &vecOuter{pos - pointingCluster.GetOuterVertex().GetPosition()};
            if (vecOuter.GetMagnitude() <= 10.f)
            {
                const float outerDot{outerDir.GetDotProduct(vecOuter)};
                if (outerDot >= -0.5f && outerDot <= 10.f)
                {
                    const float outerCross{outerDir.GetCrossProduct(vecOuter).GetMagnitude()};
                    if (outerCross < 1.f)
                    {
                        outerClusterToCandidateMap[pCluster].emplace_back(pCaloHit);
                        if (hitAssocMap.find(pCaloHit) == hitAssocMap.end())
                            hitAssocMap[pCaloHit] = 1;
                        else
                            ++hitAssocMap[pCaloHit];
                    }
                }
                //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pCaloHit->GetPositionVector(), "hit", BLACK, 1));
                //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }
        }

        CaloHitList &innerHits{innerClusterToCandidateMap[pCluster]};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &innerHits, "inner", ORANGE));
        CaloHitList &outerHits{outerClusterToCandidateMap[pCluster]};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &outerHits, "outer", GREEN));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        // Given candidates, if there are no competing pointing clusters, can probably just merge them
        // If there are competing candidates, need a more careful analysis of how hits should be allocated
        // After merging the unambiguous cases, can re-evaluate options. Consider if the competing cases
        // form crossing tracks and if so, split at the intersection, if not, perhaps longest track wins etc
        // Can use the dot product to figure out the order to walk through hits
        // Look for continuity in hits (if gaps, demand stricter perpendicular distance)
        //
        // Get CartesianPointVector, suitably sorted? and perform a sliding linear fit. Can examine the fit rms along the fit
        // and decide whether or not a hit should be added
        // Can probably collect sliding fit results and compare them to decide how to handle ambiguity
        // When stopping, probably need to get the local position from the longitudinal position and find the nearest hit
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexHitAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapLayers", m_maxGapLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FitLayers", m_fitLayers));

    float maxGapDistance = std::sqrt(m_maxGapDistanceSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapDistance", maxGapDistance));
    m_maxGapDistanceSquared = maxGapDistance * maxGapDistance;

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitSizeZ", m_hitSizeZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitSizeX", m_hitSizeX));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
