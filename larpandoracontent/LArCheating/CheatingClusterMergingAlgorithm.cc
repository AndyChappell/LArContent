/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.h"
#include "Helpers/MCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterMergingAlgorithm::CheatingClusterMergingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterMergingAlgorithm::Run()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    if (!pMCParticleList || pMCParticleList->empty())
        return STATUS_CODE_NOT_FOUND;

    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));
    if (!pClusterList || pClusterList->empty())
        return STATUS_CODE_SUCCESS;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MergeClusters(pClusterList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterMergingAlgorithm::MergeClusters(const ClusterList *pClusterList)
{
    std::map<const Cluster *, ClusterVector> matchedClusters;
    std::map<const Cluster *, bool> usedClusters;
    std::map<const Cluster *, bool> baseClusters;
    for (auto baseIter = pClusterList->begin(); baseIter != pClusterList->end(); ++baseIter)
    {
        const Cluster *pBaseCluster{*baseIter};
        if (usedClusters.find(pBaseCluster) != usedClusters.end())
            continue;
        try
        {
            // ATTN: Not all hits have populated MC weight lists, so small clusters may not have matching MC
            const MCParticle *pBaseMC(MCParticleHelper::GetMainMCParticle(pBaseCluster));
            if (!pBaseMC)
                continue;
            usedClusters[pBaseCluster] = true;
            baseClusters[pBaseCluster] = true;

            for (auto compIter = std::next(baseIter); compIter != pClusterList->end(); ++compIter)
            {
                const Cluster *pCompCluster{*compIter};
                if (usedClusters.find(pCompCluster) != usedClusters.end())
                    continue;

                try
                {
                    const MCParticle *pCompMC{MCParticleHelper::GetMainMCParticle(pCompCluster)};
                    if (!pCompMC)
                        return STATUS_CODE_FAILURE;
                    if (pBaseMC == pCompMC)
                    {
                        matchedClusters[pBaseCluster].emplace_back(pCompCluster);
                        usedClusters[pCompCluster] = true;
                    }
                }
                catch (const StatusCodeException &)
                {
                    // ATTN: While not actually used, we can avoid pointless lookups later
                    usedClusters[pCompCluster] = true;
                    continue;
                }
            }
        }
        catch (const StatusCodeException &)
        {
            continue;
        }
    }

    for (const Cluster *pCluster : *pClusterList)
    {
        if (matchedClusters.find(pCluster) != matchedClusters.end())
        {
            for (const Cluster *pMatchedCluster : matchedClusters[pCluster])
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pCluster,
                    pMatchedCluster, m_clusterListName, m_clusterListName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

