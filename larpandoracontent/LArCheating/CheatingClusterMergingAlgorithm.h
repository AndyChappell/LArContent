/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_CLUSTER_MERGING_ALGORITHM_H
#define LAR_CHEATING_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  CheatingClusterMergingAlgorithm class
 */
class CheatingClusterMergingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingClusterMergingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Merge clusters based on common best MC particle
     *
     *  @param  pClusterList the input list of clusters
     *
     *  @return The StatusCode
     */
    pandora::StatusCode MergeClusters(const pandora::ClusterList *pClusterList);


    std::string     m_clusterListName;          ///< The cluster list name
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CLUSTER_MERGING_ALGORITHM_H
