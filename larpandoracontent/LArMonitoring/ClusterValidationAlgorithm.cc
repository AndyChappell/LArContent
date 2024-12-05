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
    m_nContributions{0},
    m_wRecoHits{0},
    m_purity{0},
    m_completeness{0},
    m_fragmentationFraction{0}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ClusterValidationAlgorithm::ClusterValidationAlgorithm() :
    m_visualize{false},
    m_writeFile{false},
    m_caloHitListName{"CaloHitList2D"}
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
    this->GetMetrics(viewClustersMap, metricsMap);

    if (m_writeFile)
    {
        for (const auto &[pCluster, metrics] : metricsMap)
        {
            int view{LArClusterHelper::GetClusterHitType(pCluster)};
            const CartesianVector &mom{metrics.m_pMC->GetMomentum()};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "view", view));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "pdg", metrics.m_pMC->GetParticleId()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "tier", LArMCParticleHelper::GetHierarchyTier(metrics.m_pMC)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_x", mom.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_y", mom.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_z", mom.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_contribs", metrics.m_nContributions));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "reco_hits", metrics.m_wRecoHits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "purity", metrics.m_purity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "completeness", metrics.m_completeness));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "fragmentation", metrics.m_fragmentationFraction));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
        }
    }

    if (m_visualize)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::GetMetrics(const ViewClustersMap &viewClustersMap, ClusterValidationAlgorithm::ClusterMetricsMap &metricsMap) const
{
    const CaloHitList *pFullCaloHitList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pFullCaloHitList));
    for (const auto &[view, clusterList] : viewClustersMap)
    {
        CaloHitList viewCaloHitList;
        for (const CaloHit *const pCaloHit : *pFullCaloHitList)
            if (pCaloHit->GetHitType() == view)
                viewCaloHitList.emplace_back(pCaloHit);
        for (const Cluster *const pCluster : clusterList)
        {
            MCParticleWeightMap clusterMCWeightMap;
            const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
            CaloHitList caloHits;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
            caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
            float totalClusterWeight{0.f};
            for (const CaloHit *const pCaloHit : caloHits)
            {
                const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
                for (const auto &[pMC, weight] : weightMap)
                {
                    if (clusterMCWeightMap.find(pMC) == clusterMCWeightMap.end())
                        clusterMCWeightMap[pMC] = 0.f;
                    clusterMCWeightMap[pMC] += weight;
                    totalClusterWeight += weight;
                }
            }
            if (totalClusterWeight > 0)
                for (const auto &[pMC, weight] : clusterMCWeightMap)
                    clusterMCWeightMap[pMC] /= totalClusterWeight;
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
                float wMatchedHits{0}, wTotalRecoHits{0};
                for (const CaloHit *const pCaloHit : caloHits)
                {
                    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
                    if (weightMap.find(pMainMC) != weightMap.end())
                    {
                        wMatchedHits += weightMap.at(pMainMC);
                    }
                    if (!weightMap.empty())
                        wTotalRecoHits += 1;
                }
                if (wTotalRecoHits > 0)
                {
                    metricsMap[pCluster].m_purity = wMatchedHits / wTotalRecoHits;
                    metricsMap[pCluster].m_pMC = pMainMC;
                    metricsMap[pCluster].m_wRecoHits = wTotalRecoHits;
                }
                float wTotalMCHits{0};
                for (const CaloHit *const pCaloHit : viewCaloHitList)
                {
                    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
                    if (weightMap.find(pMainMC) != weightMap.end())
                        wTotalMCHits += weightMap.at(pMainMC);
                }

                if (wTotalMCHits > 0)
                    metricsMap[pCluster].m_completeness = wMatchedHits / wTotalMCHits;
                float accountedWeight{0.f};
                for (const auto &[pMC, weight] : clusterMCWeightMap)
                {
                    if (weight > 0.2f)
                    {
                        ++metricsMap[pCluster].m_nContributions;
                        accountedWeight += weight;
                    }
                }
                metricsMap[pCluster].m_fragmentationFraction = 1.f - accountedWeight;
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

