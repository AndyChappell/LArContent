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

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EnvelopeAssociationAlgorithm class
 */
class EnvelopeAssociationAlgorithm : public ClusterMergingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EnvelopeAssociationAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Identifies good seed clusters and then attempts to grow them based on PCA-derived bounding regions
     *
     *  @param  clusterVector the vector of clusters to consider for association
     *  @param  clusterMergeMap the map from seed clusters to candidate clusters for growth
     *  @param  minSeedCaloHits the minimum number of calo hits for a cluster to be considered a seed
     */
    void AssociateClusters(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap, const unsigned int minSeedCaloHits) const;

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
     *  @brief  Determines if a point is contained within a triangle using the Barycentric technique, as described in
     *          Real-Time Collision Detection, C. Ericson (Morgan Kaufmann, 2005).
     *
     *  @param  ab the vector describing the triangle edge ab
     *  @param  ac the vector describing the traingle edge ac
     *  @param  ap the vector describing the displacement of the point p from the point a
     *
     *  @return true if the point is contained within the triangle, false otherwise
     */
    bool IsPointContained(const pandora::CartesianVector &ab, const pandora::CartesianVector &ac, const pandora::CartesianVector &ap) const;

    /**
     *  @brief  Retrieve the vertices for the bounding envelopes to be used during cluster association
     *
     *  @param  pCluster the cluster for which the bounding envelope should be constructed
     *  @param  boundingVertices the output vector of vertices describing the bounding shapes
     */
    void GetBoundingShapes(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &boundingVertices) const;

    unsigned int m_minClusterLayers;    ///< minimum allowed number of layers for a clean cluster
    unsigned int m_minSeedCaloHits;     ///< minimum number of calo hits to form a seed for PCA
    mutable int m_runCount;             ///< the number of times the algorithm has been run
    bool m_visualize;                   ///< whether or not to visualize the algorithm
};

} // namespace lar_content

#endif // #ifndef LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H
