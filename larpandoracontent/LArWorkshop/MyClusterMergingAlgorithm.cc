/**
 *  @file   larpandoracontent/LArWorkshop/MyClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of a custom cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArWorkshop/MyClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MyClusterMergingAlgorithm::MyClusterMergingAlgorithm() : m_maxClusterSeparation{5.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusterMergingAlgorithm::Run()
{
    ClusterMergeMap clusterAssocMap;
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList));

    if (pClusterList->empty())
        return STATUS_CODE_SUCCESS;

    for (auto iter1 = pClusterList->begin(); iter1 != pClusterList->end(); ++iter1)
    {
        const Cluster *const pCluster1{*iter1};
        for (auto iter2 = std::next(iter1); iter2 != pClusterList->end(); ++iter2)
        {
            const Cluster *const pCluster2{*iter2};
            if (this->AreClustersAssociated(pCluster1, pCluster2))
            {
                clusterAssocMap[pCluster1].emplace_back(pCluster2);
                clusterAssocMap[pCluster2].emplace_back(pCluster1);
            }
        }
    }

    this->MergeClusters(clusterAssocMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MyClusterMergingAlgorithm::AreClustersAssociated(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    if (LArClusterHelper::GetClosestDistance(pCluster1, pCluster2) > m_maxClusterSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusterMergingAlgorithm::MergeClusters(const ClusterMergeMap &clusterMergeMap)
{
    ClusterSet usedClusters;
    // Need to be careful not to go into infinite loops because of bi-directional associations
    for (const auto &[pSeedCluster, associatedClusterList] : clusterMergeMap)
    {
        ClusterSet clustersToMerge;
        if (usedClusters.find(pSeedCluster) != usedClusters.end())
        {
            continue;
        }

        usedClusters.insert(pSeedCluster);
        this->TraceAssociations(clusterMergeMap, pSeedCluster, usedClusters, clustersToMerge);
        for (const Cluster *pCluster : clustersToMerge)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pCluster,
                m_inputClusterListName, m_inputClusterListName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MyClusterMergingAlgorithm::TraceAssociations(const ClusterMergeMap &clusterMergeMap, const Cluster *const pSeedCluster,
    ClusterSet &usedClusters, ClusterSet &clustersToMerge)
{
    for (const Cluster *pCluster : clusterMergeMap.at(pSeedCluster))
    {
        if (usedClusters.find(pCluster) != usedClusters.end())
        {
            continue;
        }
        else
        {
            clustersToMerge.insert(pCluster);
            usedClusters.insert(pCluster);
            this->TraceAssociations(clusterMergeMap, pCluster, usedClusters, clustersToMerge);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));
    std::cout << "Cluster list name " << m_inputClusterListName << std::endl;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterSeparation", m_maxClusterSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

