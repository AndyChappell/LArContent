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
    m_minClusterLayers{1},
    m_minSeedCaloHits{30},
    m_runCount{0},
    m_visualize{false}
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

void EnvelopeAssociationAlgorithm::PopulateClusterMergeMap(const ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    if (m_runCount > 1)
        return;
    if (m_runCount == 0)
        this->AssociateClusters(clusterVector, clusterMergeMap, m_minSeedCaloHits);
    else
        this->AssociateClusters(clusterVector, clusterMergeMap, 3 * m_minSeedCaloHits);
    ++m_runCount;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnvelopeAssociationAlgorithm::AssociateClusters(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap,
    const unsigned int minSeedCaloHits) const
{
    ClusterList seedClusters, targetClusters;
    for (ClusterVector::const_iterator iter = clusterVector.begin(); iter != clusterVector.end(); ++iter)
    {
        const Cluster *const pCluster{*iter};
        if (pCluster->GetNCaloHits() > minSeedCaloHits)
            seedClusters.emplace_back(pCluster);
        else
            targetClusters.emplace_back(pCluster);
    }
    for (const Cluster *pSeed : seedClusters)
    {
        CartesianPointVector boundingVertices;
        this->GetBoundingShapes(pSeed, boundingVertices);

        // Find the clusters with a significant fraction of hits (50%) contained by the cone and associate
        for (const Cluster * pTarget : targetClusters)
        {
            if (clusterMergeMap.find(pTarget) != clusterMergeMap.end())
                continue;
            if (!this->IsClusterContained(boundingVertices, pTarget))
                continue;

            clusterMergeMap[pSeed].emplace_back(pTarget);
            clusterMergeMap[pTarget].emplace_back(pSeed);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnvelopeAssociationAlgorithm::IsClusterContained(const CartesianPointVector &boundingVertices, const Cluster *const pCluster) const
{
    // Vertices of the extended cone
    const CartesianVector &a(boundingVertices[0]);
    const CartesianVector &b(boundingVertices[3]);
    const CartesianVector &c(boundingVertices[4]);
    CartesianVector ab(b); ab -= a;
    CartesianVector ac(c); ac -= a;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);
    CaloHitList containedHits, uncontainedHits;
    for (const CaloHit *pCaloHit : caloHits)
    {
        CartesianVector ap(pCaloHit->GetPositionVector()); ap -= a;
        if (this->IsPointContained(ab, ac, ap))
            containedHits.emplace_back(pCaloHit);
        else
            uncontainedHits.emplace_back(pCaloHit);
    }

    return !containedHits.empty() && containedHits.size() >= uncontainedHits.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnvelopeAssociationAlgorithm::IsPointContained(const CartesianVector &ab, const CartesianVector &ac, const CartesianVector &ap) const
{
    const float dotABAB{ab.GetDotProduct(ab)};
    const float dotABAP{ab.GetDotProduct(ap)};
    const float dotACAC{ac.GetDotProduct(ac)};
    const float dotACAB{ac.GetDotProduct(ab)};
    const float dotACAP{ac.GetDotProduct(ap)};

    const float denom{dotACAC * dotABAB - dotACAB * dotACAB};
    if (std::fabs(denom) < std::numeric_limits<float>::epsilon())
        return false;
    const float invDenom{1.f / denom};
    const float u{(dotABAB * dotACAP - dotACAB * dotABAP) * invDenom};
    const float v{(dotACAC * dotABAP - dotACAB * dotACAP) * invDenom};

    return (u >= 0) && (v >= 0) && ((u + v) <= 1);
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

    CartesianVector p2(p1);
    CartesianVector p2r(eigenVectors[1]), p2l(0.f, 0.f, 0.f);
    CartesianVector p3r(eigenVectors[1]), p3l(0.f, 0.f, 0.f);
    if (m_runCount > 0)
    {   // Broader cone surrounding the seed cluster
        p2r *= 3.5f * length * eigenRatio;
        p2l -= p2r;
        p2r += p2; p2l += p2;
        p3r = p2r; p3l = p2l;
    }
    else
    {
        // Define an extended cone beyond the end of the seed cluster
        p2 += dl;
        p2r *= 2 * length * eigenRatio;
        p2l -= p2r;
        p2r += p2; p2l += p2;

        // Define an interior box within the extended cone
        p3r *= length * eigenRatio;
        p3l -= p3r;
        p3r += p2; p3l += p2;
    }

    if (m_visualize && LArClusterHelper::GetClusterHitType(pCluster) == HitType::TPC_VIEW_W)
    {
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p1r, "A", BLUE, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p1l, "B", BLUE, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p0, "C", BLUE, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p2r, "D", RED, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p2l, "E", RED, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1r, &p3r, "F", GREEN, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p1l, &p3l, "G", GREEN, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p2r, &p2l, "H", RED, 3, 1));
        PANDORA_MONITORING_API(Pause(this->GetPandora()));
    }

    if (m_visualize && LArClusterHelper::GetClusterHitType(pCluster) == HitType::TPC_VIEW_W)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, "Seed", BLUE));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p2l, "AB", BLACK, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p2r, "AC", BLACK, 3, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p2l, &p2r, "BC", BLACK, 3, 1));
        PANDORA_MONITORING_API(Pause(this->GetPandora()));
    }

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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers",
        m_minClusterLayers));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSeedCaloHits",
        m_minSeedCaloHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize",
        m_visualize));

    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
