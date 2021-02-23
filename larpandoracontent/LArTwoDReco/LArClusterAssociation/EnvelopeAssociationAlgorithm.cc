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
    for (ClusterVector::const_iterator iter1 = clusterVector.begin(); iter1 != clusterVector.end(); ++iter1)
    {
        const Cluster *const pCluster1{*iter1};
        // Remove hard-coding
        if (pCluster1->GetNCaloHits() > 30)
        {
            // Get the eigen vectors for this cluster
            const OrderedCaloHitList &orderedCaloHitList(pCluster1->GetOrderedCaloHitList());
            CaloHitList caloHits;
            orderedCaloHitList.FillCaloHitList(caloHits);
            CartesianVector centroid(0.f, 0.f, 0.f);
            LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVectors;
            LArPcaHelper::RunPca(caloHits, centroid, eigenValues, eigenVectors);

            // Construct a cone around the hits based on the eigen vectors
            CartesianVector p0(0.f, 0.f, 0.f), p1(0.f, 0.f, 0.f);
            LArClusterHelper::GetExtremalCoordinates(pCluster1, p0, p1);
            p0 -= centroid;
            p1 -= centroid;
            const float p0Dote0Norm{p0.GetDotProduct(eigenVectors[0]) / eigenVectors[0].GetMagnitudeSquared()};
            const float p1Dote0Norm{p1.GetDotProduct(eigenVectors[0]) / eigenVectors[0].GetMagnitudeSquared()};
            p0.SetValues(p0Dote0Norm * eigenVectors[0].GetX(), p0Dote0Norm * eigenVectors[0].GetY(), p0Dote0Norm * eigenVectors[0].GetZ());
            p0 += centroid;
            p1.SetValues(p1Dote0Norm * eigenVectors[0].GetX(), p1Dote0Norm * eigenVectors[0].GetY(), p1Dote0Norm * eigenVectors[0].GetZ());
            p1 += centroid;
            CartesianVector dl(p1);
            dl -= p0;
            const float length{dl.GetMagnitude()};
            CartesianVector b(eigenVectors[1]), c(0.f, 0.f, 0.f);
            b *= length * std::sqrt(eigenValues.GetY() / eigenValues.GetX());
            c -= b;
            b += p1; c+= p1;
            CartesianVector p2(p1);
            p2 += dl;
            CartesianVector d(eigenVectors[1]), e(0.f, 0.f, 0.f);
            d *= 2 * length * std::sqrt(eigenValues.GetY() / eigenValues.GetX());
            e -= d;
            d += p2; e += p2;
            CartesianVector f(eigenVectors[1]), g(0.f, 0.f, 0.f);
            f *= length * std::sqrt(eigenValues.GetY() / eigenValues.GetX());
            g -= f;
            f += p2; g += p2; 

            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &b, "A", BLUE, 3, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &b, &c, "B", BLUE, 3, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &c, &p0, "C", BLUE, 3, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &b, &d, "D", BLUE, 3, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &c, &e, "E", BLUE, 3, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &b, &f, "F", BLUE, 3, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &c, &g, "G", BLUE, 3, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &d, &e, "H", BLUE, 3, 1));
            PANDORA_MONITORING_API(Pause(this->GetPandora()));
        }

        // Find the clusters with a significant fraction of hits (50%) contained by the cone and associate
        for (ClusterVector::const_iterator iter2 = std::next(iter1); iter2 != clusterVector.end(); ++iter2)
        {
            const Cluster *const pCluster2{*iter2};

            if (!this->AreClustersAssociated(pCluster1, pCluster2))
                continue;

            clusterAssociationMap[pCluster1].m_forwardAssociations.insert(pCluster2);
            clusterAssociationMap[pCluster2].m_backwardAssociations.insert(pCluster1);
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

bool EnvelopeAssociationAlgorithm::AreClustersAssociated(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    // This is inefficient, considering computing all centroids first and looking them up
    CartesianVector centroid1(0.f, 0.f, 0.f);
    if (LArClusterHelper::GetCentroid(pCluster1, centroid1) != STATUS_CODE_SUCCESS)
        return false;
    CartesianVector centroid2(0.f, 0.f, 0.f);
    if (LArClusterHelper::GetCentroid(pCluster2, centroid2) != STATUS_CODE_SUCCESS)
        return false;

    return centroid1.GetDistanceSquared(centroid2) <= m_maxGapDistanceSquared;
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
