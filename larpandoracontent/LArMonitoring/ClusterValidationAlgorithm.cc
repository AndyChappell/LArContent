/**
 *  @file   larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

ClusterValidationAlgorithm::ClusterMetrics::ClusterMetrics() :
    m_pMC{nullptr},
    m_nRecoHits{0},
    m_purity{0},
    m_completeness{0}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ClusterValidationAlgorithm::ClusterValidationAlgorithm() :
    m_visualize{false},
    m_writeFile{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterValidationAlgorithm::~ClusterValidationAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterValidationAlgorithm::Run()
{
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    }

    ViewClustersMap viewClustersMap;
    for (std::string listName : m_clusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pClusterList));
        for (const Cluster *const pCluster : *pClusterList)
        {
            const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
            viewClustersMap[view].emplace_back(pCluster);
        }
    }

    ClusterMetricsMap metricsMap;
    this->GetPurity(viewClustersMap, metricsMap);

    for (const auto &[pCluster, metrics] : metricsMap)
    {
        std::cout << "Cluster: " << pCluster << " nHits: " << metrics.m_nRecoHits << " purity: " << metrics.m_purity << std::endl;
    }

    if (m_visualize)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::GetPurity(const ViewClustersMap &viewClustersMap, ClusterValidationAlgorithm::ClusterMetricsMap &metricsMap) const
{
    std::cout << "Map size: " << viewClustersMap.size() << std::endl;
    for (const auto &[view, clusterList] : viewClustersMap)
    {
        std::cout << "Num clusters: " << clusterList.size() << std::endl;
        for (const Cluster *const pCluster : clusterList)
        {
            MCParticleWeightMap clusterMCWeightMap;
            const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
            CaloHitList caloHits;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
            caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
            for (const CaloHit *const pCaloHit : caloHits)
            {
                const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
                for (const auto &[pMC, weight] : weightMap)
                {
                    if (clusterMCWeightMap.find(pMC) == clusterMCWeightMap.end())
                        clusterMCWeightMap[pMC] = 0.f;
                    clusterMCWeightMap[pMC] += weight;
                }
            }
            // Determine which MC particle contributes the most weight across the cluster
            const MCParticle *pMainMC{metricsMap.find(pCluster) != metricsMap.end() ? metricsMap[pCluster].m_pMC : nullptr};
            if (!pMainMC)
            {
                float maxWeight{std::numeric_limits<float>::lowest()};
                for (const auto &[pMC, weight] : clusterMCWeightMap)
                {
                    if (weight > maxWeight)
                    {
                        pMainMC = pMC;
                        maxWeight = weight;
                    }
                }
            }
            if (pMainMC)
            {
                int nMatchedHits{0}, nTotalHits{0};
                for (const CaloHit *const pCaloHit : caloHits)
                {
                    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
                    if (weightMap.find(pMainMC) != weightMap.end())
                        ++nMatchedHits;
                    if (!weightMap.empty())
                        ++nTotalHits;
                }
                if (nTotalHits)
                {
                    metricsMap[pCluster].m_purity = nMatchedHits / static_cast<float>(nTotalHits);
                    metricsMap[pCluster].m_pMC = pMainMC;
                    metricsMap[pCluster].m_nRecoHits = nTotalHits;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
