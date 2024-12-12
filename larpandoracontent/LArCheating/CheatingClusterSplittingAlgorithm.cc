/**
 *  @file   larpandoracontent/LArMonitoring/CheatingClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingClusterSplittingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterSplittingAlgorithm::CheatingClusterSplittingAlgorithm() :
    m_visualize{false},
    m_caloHitListName{"CaloHitList2D"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CheatingClusterSplittingAlgorithm::~CheatingClusterSplittingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterSplittingAlgorithm::Run()
{
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    }
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList2D));
    MCParticleWeightMap allMCMap;
    for (const CaloHit *const pCaloHit : *pCaloHitList2D)
    {
        const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
        for (const auto &[pMC, weight] : weightMap)
        {
            if (allMCMap.find(pMC) == allMCMap.end())
                allMCMap[pMC] = 0.f;
            allMCMap[pMC] += weight;
        }
    }

    std::map<const Cluster *, std::set<CaloHitList>> clusterSplitMap;
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));
    for (const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
        const CaloHitList &isolatedHitList{pCluster->GetIsolatedCaloHitList()};
        caloHitList.insert(caloHitList.end(), isolatedHitList.begin(), isolatedHitList.end());
        MCParticleWeightMap clusterMCMap;
        float totalClusterWeight{0.f};
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
            for (const auto &[pMC, weight] : weightMap)
            {
                if (clusterMCMap.find(pMC) == clusterMCMap.end())
                    clusterMCMap[pMC] = 0.f;
                clusterMCMap[pMC] += weight;
                totalClusterWeight += weight;
            }
        }
        if (totalClusterWeight > 0)
        {
            for (const auto &[pMC, weight] : clusterMCMap)
                clusterMCMap[pMC] /= totalClusterWeight;

            float maxWeight{0.f};
            const MCParticle *pMainMC{nullptr};
            std::set<const MCParticle *> mcSet;
            for (const auto &[pMC, weight] : clusterMCMap)
            {
                if ((weight > 0.2f) || (allMCMap.at(pMC) > 0.2))
                    mcSet.insert(pMC);
                if (weight > maxWeight)
                {
                    pMainMC = pMC;
                    maxWeight = weight;
                }
            }
            if (mcSet.size() >  1)
            {
                mcSet.erase(pMainMC);
                LArMCParticleHelper::MCContributionMap mcHitMap;
                for (const CaloHit *pCaloHit : caloHitList)
                {
                    try
                    {
                        const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                        if (mcSet.find(pMC) != mcSet.end())
                            mcHitMap[pMC].emplace_back(pCaloHit);
                        else
                            mcHitMap[pMainMC].emplace_back(pCaloHit);
                    }
                    catch(StatusCodeException e)
                    {
                        // Not a real MC mapping, but we just want to ensure any hits without MC end up in the "original" cluster
                        mcHitMap[pMainMC].emplace_back(pCaloHit);
                    }
                }
                for (const auto &[pMC, mcCaloHitList] : mcHitMap)
                    clusterSplitMap[pCluster].insert(mcCaloHitList);
            }
        }
    }

    for (const auto &[pCluster, clusterHitsSet] : clusterSplitMap)
    {
        const ClusterList clusterList(1, pCluster);
        std::string clusterListToSaveName, clusterListToDeleteName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList,
            clusterListToDeleteName, clusterListToSaveName));
        for (const CaloHitList &caloHitList : clusterHitsSet)
        {
            const Cluster *pNewCluster(nullptr);
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList = caloHitList;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
        }
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
    }

    if (m_visualize)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

