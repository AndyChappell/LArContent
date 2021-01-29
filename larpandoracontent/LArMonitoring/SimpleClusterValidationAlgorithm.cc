/**
 *  @file   larpandoracontent/LArMonitoring/SimpleClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the simple cluster validation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/SimpleClusterValidationAlgorithm.h"
#include "Helpers/MCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

SimpleClusterValidationAlgorithm::SimpleClusterValidationAlgorithm() :
    m_writeFile{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SimpleClusterValidationAlgorithm::~SimpleClusterValidationAlgorithm()
{
    if (m_writeFile)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_filename, "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "SimpleClusterValidationAlgorithm: Unable to write tree to file " << m_filename << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterValidationAlgorithm::Run()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    if (!pMCParticleList || pMCParticleList->empty())
        return STATUS_CODE_NOT_FOUND;

    const ClusterList *pClusterList{nullptr};
    for (std::string clusterListName : m_clusterListNames)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));
        if (!pClusterList || pClusterList->empty())
            continue;

        LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
        LArMCParticleHelper::MCContributionMap mcToHitMap;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetRecoTruthMaps(pMCParticleList, pClusterList, hitToMCMap, mcToHitMap));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GenerateMetrics(pClusterList, hitToMCMap, mcToHitMap));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterValidationAlgorithm::GenerateMetrics(const ClusterList *pClusterList,
    const LArMCParticleHelper::CaloHitToMCMap &hitToMCMap, const LArMCParticleHelper::MCContributionMap &mcToHitMap)
{
    for (const Cluster *pCluster : *pClusterList)
    {
        const MCParticle *pMC{nullptr};
        try
        {
            // ATTN: Not all hits have a populated MC particle weight list, so small clusters will sometimes not have an MC match
            pMC = MCParticleHelper::GetMainMCParticle(pCluster);
            if (!pMC)
                continue;
        }
        catch(const StatusCodeException &)
        {
            continue;
        }
        if (mcToHitMap.find(pMC) == mcToHitMap.end())
        {
            std::cout << "Found an MC particle not in the MC to hits map" << std::endl;
            continue;
        }
        const CaloHitList &mcHits{mcToHitMap.at(pMC)};
        const int nMCHits{static_cast<int>(mcHits.size())};

        CaloHitList clusterHits;
        for (const auto &clusterHitPair : pCluster->GetOrderedCaloHitList())
        {
            CaloHitList caloHits(*clusterHitPair.second);
            clusterHits.merge(caloHits);
        }
        CaloHitList isoHits(pCluster->GetIsolatedCaloHitList());
        clusterHits.merge(isoHits);

        const int view{static_cast<int>(clusterHits.front()->GetHitType())};
        const int nClusterHits{static_cast<int>(clusterHits.size())};
        int nSharedHits{0};
        int pdg{0};
        int isTrueTrack(1);
        for (const CaloHit *pCaloHit : clusterHits)
        {
            const auto iter(hitToMCMap.find(pCaloHit));
            if (iter != hitToMCMap.end())
            {
                const MCParticle *pMatchedMC{iter->second};
                if (pMatchedMC == pMC)
                {
                    ++nSharedHits;
                    pdg = pMC->GetParticleId();
                    if (std::abs(pdg) == PHOTON || std::abs(pdg) == E_MINUS)
                        isTrueTrack = 0;
                }
            }
        }
        const float purity{static_cast<float>(nSharedHits) / nClusterHits};
        const float completeness{static_cast<float>(nSharedHits) / nMCHits};

        if (m_writeFile)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "view", view));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "nClusterHits", nClusterHits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "nMCHits", nMCHits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "nSharedHits", nSharedHits));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "purity", purity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "completeness", completeness));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "pdg", pdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "isTrueTrack", isTrueTrack));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterValidationAlgorithm::GetRecoTruthMaps(const MCParticleList *pMCList, const ClusterList *pClusterList,
    LArMCParticleHelper::CaloHitToMCMap &hitToMCMap, LArMCParticleHelper::MCContributionMap &mcToHitMap)
{
    LArMCParticleHelper::MCRelationMap mcToTargetMap;
    LArMCParticleHelper::GetMCToSelfMap(pMCList, mcToTargetMap);

    CaloHitList caloHits;
    for (const Cluster *pCluster : *pClusterList)
    {
        for (const auto &clusterHitPair : pCluster->GetOrderedCaloHitList())
        {
            CaloHitList clusterHits(*clusterHitPair.second);
            caloHits.merge(clusterHits);
        }
        CaloHitList clusterIsoHits(pCluster->GetIsolatedCaloHitList());
        caloHits.merge(clusterIsoHits);
    }

    try
    {
        LArMCParticleHelper::GetMCParticleToCaloHitMatches(&caloHits, mcToTargetMap, hitToMCMap, mcToHitMap);
    }
    catch (const StatusCodeException &e)
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
        if (m_filename.empty())
            return STATUS_CODE_INVALID_PARAMETER;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
        if (m_treeName.empty())
            return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

