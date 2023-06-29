/**
 *  @file   larpandoracontent/LArWorkshop/MyClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of a custom cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArWorkshop/MyClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MyClusterMergingAlgorithm::MyClusterMergingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusterMergingAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MyClusterMergingAlgorithm::AreClustersAssociated(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    (void)pCluster1;
    (void)pCluster2;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusterMergingAlgorithm::MergeClusters(const ClusterMergeMap &clusterMergeMap) const
{
    (void)clusterMergeMap;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    (void)xmlHandle;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

