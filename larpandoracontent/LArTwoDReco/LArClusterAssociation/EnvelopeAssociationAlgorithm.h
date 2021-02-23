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
     *  @brief  Determine whether two clusters are associated
     *
     *  @param  pCluster1 the inner cluster
     *  @param  pCluster2 the outer cluster
     *
     *  @return whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    unsigned int m_minClusterLayers;            ///< minimum allowed number of layers for a clean cluster
    float        m_maxGapDistanceSquared;       ///< maximum allowed distance (squared) between associated clusters
};

} // namespace lar_content

#endif // #ifndef LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H
