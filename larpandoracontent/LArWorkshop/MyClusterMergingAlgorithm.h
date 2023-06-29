/**
 *  @file   larpandoracontent/LArWorkshop/MyClusterMergingAlgorithm.h
 *
 *  @brief  Header file for a custom cluster merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MY_CLUSTER_MERGING_ALGORITHM_H
#define LAR_MY_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MyClusterMergingAlgorithm class
 */
class MyClusterMergingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MyClusterMergingAlgorithm();

private:
    typedef std::map<const pandora::Cluster *, pandora::ClusterList> ClusterMergeMap;

    pandora::StatusCode Run();

    /**
     *  @brief Determine if two clusters are associated with each other
     *
     *  @param pCluster1 The first cluster to consider in the candidate association
     *  @param pCluster2 The second cluster to consider in the candidate association
     *
     *  @return Whether or not the pair of clusters is associated
     */
    bool AreClustersAssociated(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief Merge together associated clusters
     *
     *  @param clusterMergeMap The map from seed ckuster to asssociated clusters describing the merges to make
     *
     *  @return The status code for the merge operation
     */
    pandora::StatusCode MergeClusters(const ClusterMergeMap &clusterMergeMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_MY_CLUSTER_MERGING_ALGORITHM_H

