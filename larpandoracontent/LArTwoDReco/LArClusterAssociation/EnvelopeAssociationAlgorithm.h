/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/EnvelopeAssociationAlgorithm.h
 *
 *  @brief  Header file for the proximity association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H
#define LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EnvelopeAssociationAlgorithm class
 */
class EnvelopeAssociationAlgorithm : public ClusterAssociationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EnvelopeAssociationAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  Determine whether a target cluster is contained by a bounding region
     *
     *  @param  boundVertices the bounding vertices
     *  @param  pCluster the target cluster
     *
     *  @return whether the cluster is contained within the bounding region
     */
    bool IsClusterContained(const pandora::CartesianPointVector &boundingVertices, const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Retrieve the vertices for the bounding envelopes to be used during cluster association
     *
     *  @param  pCluster the cluster for which the bounding envelope should be constructed
     *  @param  boundingVertices the output vector of vertices describing the bounding shapes
     */
    void GetBoundingShapes(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &boundingVertices) const;

    unsigned int m_minClusterLayers;            ///< minimum allowed number of layers for a clean cluster
    float        m_maxGapDistanceSquared;       ///< maximum allowed distance (squared) between associated clusters
};

} // namespace lar_content

#endif // #ifndef LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H
