/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/EnvelopeAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the proximity association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/EnvelopeAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EnvelopeAssociationAlgorithm::EnvelopeAssociationAlgorithm() :
    m_minClusterLayers(1),
    m_maxGapDistanceSquared(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnvelopeAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
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

void EnvelopeAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterList seedClusters, targetClusters;
    for (ClusterVector::const_iterator iter = clusterVector.begin(); iter != clusterVector.end(); ++iter)
    {
        const Cluster *const pCluster{*iter};
        // Remove hard-coding
        if (pCluster->GetNCaloHits() > 30)
            seedClusters.emplace_back(pCluster);
        else
            targetClusters.emplace_back(pCluster);
    }
    for (const Cluster *pSeed : seedClusters)
    {
        std::cout << "NHits " << pSeed << " " << pSeed->GetNCaloHits() << std::endl;
        CartesianPointVector boundingVertices;
        this->GetBoundingShapes(pSeed, boundingVertices);

        // Find the clusters with a significant fraction of hits (50%) contained by the cone and associate
        for (const Cluster * pTarget : targetClusters)
        {
            if (!this->IsClusterContained(boundingVertices, pTarget))
                continue;

            clusterAssociationMap[pSeed].m_forwardAssociations.insert(pTarget);
            clusterAssociationMap[pTarget].m_backwardAssociations.insert(pSeed);
        }
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnvelopeAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnvelopeAssociationAlgorithm::IsClusterContained(const CartesianPointVector &boundingVertices, const Cluster *const pCluster) const
{
    // Vertices of the extended cone
    const CartesianVector &a(boundingVertices[0]);
    const CartesianVector &b(boundingVertices[3]);
    const CartesianVector &c(boundingVertices[4]);

    (void)a;
    (void)b;
    (void)c;
    (void)pCluster;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnvelopeAssociationAlgorithm::GetBoundingShapes(const Cluster *const pCluster, CartesianPointVector &boundingVertices) const
{
    // Retrieve the hits for PCA
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);

    // Get the eigen vectors for this cluster
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(caloHits, centroid, eigenValues, eigenVectors);

    // Project the extremal hits of the cluster onto the principal axis
    CartesianVector p0(0.f, 0.f, 0.f), p1(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pCluster, p0, p1);
    p0 -= centroid;
    p1 -= centroid;
    const CartesianVector &principalAxis(eigenVectors[0]);
    const float vMagSquared{principalAxis.GetMagnitudeSquared()};
    const float normDot0p{p0.GetDotProduct(principalAxis) / vMagSquared};
    const float normDot1p{p1.GetDotProduct(principalAxis) / vMagSquared};
    p0.SetValues(normDot0p * principalAxis.GetX() + centroid.GetX(), normDot0p * principalAxis.GetY() + centroid.GetY(),
        normDot0p * principalAxis.GetZ() + centroid.GetZ());
    p1.SetValues(normDot1p * principalAxis.GetX() + centroid.GetX(), normDot1p * principalAxis.GetY() + centroid.GetY(),
        normDot1p * principalAxis.GetZ() + centroid.GetZ());

    // Get the length of the cluster along the principal axis
    CartesianVector dl(p1); dl -= p0;
    const float length{dl.GetMagnitude()};

    const float eigenRatio{std::sqrt(eigenValues.GetY() / eigenValues.GetX())};
    // Define a cone that (approximately) envelopes the cluster
    CartesianVector p1r(eigenVectors[1]), p1l(0.f, 0.f, 0.f);
    p1r *= length * eigenRatio;
    p1l -= p1r;
    p1r += p1; p1l += p1;

    // Define an extended cone beyond the end of the cluster
    CartesianVector p2(p1); p2 += dl;
    CartesianVector p2r(eigenVectors[1]), p2l(0.f, 0.f, 0.f);
    p2r *= 2 * length * eigenRatio;
    p2l -= p2r;
    p2r += p2; p2l += p2;

    // Define an interior box within the extended cone
    CartesianVector p3r(eigenVectors[1]), p3l(0.f, 0.f, 0.f);
    p3r *= length * eigenRatio;
    p3l -= p3r;
    p3r += p2; p3l += p2; 

    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p1r, "A", BLUE, 3, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p1l, "B", BLUE, 3, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p0, "C", BLUE, 3, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p2r, "D", RED, 3, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p2l, "E", RED, 3, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p3r, "F", GREEN, 3, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p3l, "G", GREEN, 3, 1));
    PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p2r, &p2l, "H", RED, 3, 1));
    PANDORA_MONITORING_API(Pause(this->GetPandora()));

    boundingVertices.emplace_back(p0);     // Origin of cone
    boundingVertices.emplace_back(p1l);    // 'Left' vertex of primary
    boundingVertices.emplace_back(p1r);    // 'Right' vertex of primary
    boundingVertices.emplace_back(p2l);    // 'Left' vertex of extended cone
    boundingVertices.emplace_back(p2r);    // 'Right' vertex of extended cone
    boundingVertices.emplace_back(p3l);    // 'Left' downstream vertex of extended box
    boundingVertices.emplace_back(p3r);    // 'Right' downstream vertex of extended box
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnvelopeAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    float maxGapDistance{std::sqrt(m_maxGapDistanceSquared)};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGapDistance", maxGapDistance));
    m_maxGapDistanceSquared = maxGapDistance * maxGapDistance;

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
