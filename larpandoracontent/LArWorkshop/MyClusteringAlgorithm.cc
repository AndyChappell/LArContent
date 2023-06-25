/**
 *  @file   larpandoracontent/LArWorkshop/MyClusteringAlgorithm.cc
 *
 *  @brief  Implementation of a custom clustering algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArWorkshop/MyClusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MyClusteringAlgorithm::MyClusteringAlgorithm() : m_clusterSize{10}, m_outputListPrefix{"Clusters"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusteringAlgorithm::Run()
{
    // lambda function to sort the CaloHit/distance pairs in ascending order of distance
    auto sort_func{[](const std::pair<const CaloHit *, float> &pair1, const std::pair<const CaloHit *, float> &pair2) -> bool
    {
        return pair1.second < pair2.second;
    }};

    // Loop over each CaloHitList
    for (std::string listName : m_caloHitListNames)
    {
        // Make a temporary list to store our output clusters and set it as the current target list
        const ClusterList *pTemporaryList{nullptr};
        std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTemporaryList, temporaryListName));

        std::map<const CaloHit *, bool> usedCaloHitMap;
        std::vector<std::pair<const CaloHit *, float>> hitsToCluster;
        // Load the current named CaloHitList into pCaloHitList
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        if (pCaloHitList->empty())
            continue;

        // Iterate over all hits in the current hit list
        for (auto iter1 = pCaloHitList->begin(); iter1 != pCaloHitList->end(); )
        {
            const CaloHit *const pCaloHit1{*iter1};
            // Skip the hit if we've already clustered it
            if (usedCaloHitMap.find(pCaloHit1) != usedCaloHitMap.end())
            {
                ++iter1;
                continue;
            }
            // Tag the seed hit of the cluster as in use
            usedCaloHitMap[pCaloHit1] = true;
            // Iterate over all subsequent hits in the hit list - inefficient solution, but we don't care for this demo
            for (auto iter2 = std::next(iter1, 1); iter2 != pCaloHitList->end(); ++iter2)
            {
                const CaloHit *const pCaloHit2{*iter2};
                // Skip the hit is we've already clustered it
                if (usedCaloHitMap.find(pCaloHit2) != usedCaloHitMap.end())
                    continue;
                // Get the distance between this and the seed hit and add it to the current cluster if the cluster is less than the maximum
                // size, or if it is closer than an existing hit in the cluster (replacing the most distant hit)
                float distance{pCaloHit2->GetPositionVector().GetDistanceSquared(pCaloHit1->GetPositionVector())};
                const size_t clusterSize{hitsToCluster.size()};
                if (clusterSize < (m_clusterSize - 1))
                {
                    hitsToCluster.emplace_back(std::make_pair(pCaloHit2, distance));
                    if (clusterSize == (m_clusterSize - 1))
                        std::sort(hitsToCluster.begin(), hitsToCluster.end(), sort_func);
                }
                else
                {
                    if (distance < hitsToCluster.back().second)
                    {
                        hitsToCluster[m_clusterSize - 2] = std::make_pair(pCaloHit2, distance);
                        std::sort(hitsToCluster.begin(), hitsToCluster.end(), sort_func);
                    }
                }
            }
            ++iter1;

            // Create a cluster from the seed hit
            const Cluster *pCluster{nullptr};
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.emplace_back(pCaloHit1);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
            // Add the associated hits to the created cluster
            for (const auto &[pAssociatedCaloHit, distance] : hitsToCluster)
            {
                usedCaloHitMap[pAssociatedCaloHit] = true;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pAssociatedCaloHit));
            }
            // Empty out temporary vector of hits/distances
            hitsToCluster.clear();
        }
        if (!pTemporaryList->empty())
        {
            // Get the view associated with the current hit list, adapt the output list name and save the list
            std::string clusterListName{m_outputListPrefix};
            const HitType view{pCaloHitList->front()->GetHitType()};
            switch (view)
            {
                case TPC_VIEW_U:
                    clusterListName += "U";
                    break;
                case TPC_VIEW_V:
                    clusterListName += "V";
                    break;
                default:
                    clusterListName += "W";
                    break;
            }
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, clusterListName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterSize", m_clusterSize));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputListPrefix", m_outputListPrefix));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

