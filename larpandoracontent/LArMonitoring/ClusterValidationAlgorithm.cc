/**
 *  @file   larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster validation class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterValidationAlgorithm::ClusterValidationAlgorithm() :
    m_writeToTree(false),
    m_eventNumber(-1),
    m_criteria(LArMCParticleHelper::IsBeamNeutrinoFinalState)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterValidationAlgorithm::~ClusterValidationAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            std::string clusterTreeName{m_treeName + "_cluster"};
            std::string mcTreeName{m_treeName + "_mc"};
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), mcTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), clusterTreeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "ClusterValidationAlgorithm: Unable to write trees to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterValidationAlgorithm::Run()
{
    ++m_eventNumber;

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap targetMCToHitsMap;
    LArMCParticleHelper::MCContributionMap allMCToHitsMap;
    this->GetMCToHitsMaps(pMCParticleList, pCaloHitList, targetMCToHitsMap, allMCToHitsMap);

    LArMCParticleHelper::ClusterContributionMap clusterToHitsMap;
    for (std::string clusterListName : m_clusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PandoraContentApi::GetList(*this, clusterListName, pClusterList);

        if(!pClusterList)
            continue;

        LArMCParticleHelper::GetClusterToReconstructable2DHitsMap(*pClusterList, allMCToHitsMap, clusterToHitsMap);
    }
    // ATTN : Ensure all mc primaries have an entry in mcToClusterHitSharingMap, even if no associated clusters.
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({allMCToHitsMap}, mcPrimaryVector);

    this->ProcessOutput(targetMCToHitsMap, clusterToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::GetMCToHitsMaps(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList,
    LArMCParticleHelper::MCContributionMap &targetMCToHitsMap, LArMCParticleHelper::MCContributionMap &allMCToHitsMap)
{
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::PrimaryParameters parameters;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, m_criteria, targetMCToHitsMap);

        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, m_criteria, allMCToHitsMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::ProcessOutput(const LArMCParticleHelper::MCContributionMap &mcToHitsMap,
    const LArMCParticleHelper::ClusterContributionMap &clusterToHitsMap)
{
    std::string clusterTreeName{m_treeName + "_cluster"};
    std::string mcTreeName{m_treeName + "_mc"};
    MCParticleVector mcVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcToHitsMap}, mcVector);
    ClusterVector clusterVector;
    LArMonitoringHelper::GetOrderedClusterVector(clusterToHitsMap, clusterVector);

    int clusterIndex{0}, mcIndex{0};
    ClusterToIdMap clusterToIdMap;
    for (const Cluster *const pCluster : clusterVector)
        clusterToIdMap.insert(ClusterToIdMap::value_type(pCluster, ++clusterIndex));

    for (const MCParticle *const pMC : mcVector)
    {
        bool hasMatch(false);
        for (const auto [ pCluster, dummy ] : clusterToHitsMap)
        {
            CaloHitList sharedHits(LArMCParticleHelper::GetSharedHits(clusterToHitsMap.at(pCluster), mcToHitsMap.at(pMC)));
            if (!sharedHits.empty())
            {
                hasMatch = true;
                break;
            }
        }

        const CaloHitList &mcHitList{mcToHitsMap.at(pMC)};
        const int pdg{pMC->GetParticleId()};
        const int nMCHitsU{static_cast<int>(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcHitList))};
        const int nMCHitsV{static_cast<int>(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcHitList))};
        const int nMCHitsW{static_cast<int>(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcHitList))};

        if (hasMatch)
        {   // Potential match to at least 1 cluster, get details
            for (const auto [ pCluster, dummy ] : clusterToHitsMap)
            {
                CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(clusterToHitsMap.at(pCluster), mcToHitsMap.at(pMC)));
                if (sharedHitList.empty())
                    continue;

                const int clusterId(clusterToIdMap.at(pCluster));
                const CaloHitList &clusterHitList(clusterToHitsMap.at(pCluster));
                const HitType view{clusterHitList.front()->GetHitType()};
                const int nSharedHits{static_cast<int>(sharedHitList.size())};
                const int nClusterHits{static_cast<int>(clusterHitList.size())};
                const int nMCHits{static_cast<int>(LArMonitoringHelper::CountHitsByType(view, mcHitList))};

                int nMCSharedHitsU{0}, nMCSharedHitsV{0}, nMCSharedHitsW{0};
                int nContributingHitsU{0}, nContributingHitsV{0}, nContributingHitsW{0};
                if (view == TPC_VIEW_U)
                {
                    nMCSharedHitsU = nSharedHits;
                    nContributingHitsU = nClusterHits;
                }
                if (view == TPC_VIEW_V)
                {
                    nMCSharedHitsV = nSharedHits;
                    nContributingHitsV = nClusterHits;
                }
                if (view == TPC_VIEW_W)
                {
                    nMCSharedHitsW = nSharedHits;
                    nContributingHitsW = nClusterHits;
                }

                // Record cluster-centric information if it was the best match - CONSIDER: Might want to retain non-matching too
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), clusterTreeName.c_str(), "eventNumber", m_eventNumber));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), clusterTreeName.c_str(), "clusterId", clusterId));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), clusterTreeName.c_str(), "mcId", mcIndex));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), clusterTreeName.c_str(), "pdg", pdg));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), clusterTreeName.c_str(), "nMCHits", nMCHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), clusterTreeName.c_str(), "nSharedHits", nSharedHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), clusterTreeName.c_str(), "nClusterHits", nClusterHits));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), clusterTreeName.c_str()));

                // Record the MC information relating to this cluster, want all clusters, not just the best
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "eventNumber", m_eventNumber));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "clusterId", clusterId));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "mcId", mcIndex));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "pdg", pdg));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nMCHitsU", nMCHitsU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nMCHitsV", nMCHitsV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nMCHitsW", nMCHitsW));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nSharedHitsU", nMCSharedHitsU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nSharedHitsV", nMCSharedHitsV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nSharedHitsW", nMCSharedHitsW));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nContributingHitsU", nContributingHitsU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nContributingHitsV", nContributingHitsV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nContributingHitsW", nContributingHitsW));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), mcTreeName.c_str()));
            }
        }
        else
        {   // No matching clusters
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "eventNumber", m_eventNumber));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "clusterId", -1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "mcId", mcIndex));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "pdg", pdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nMCHitsU", nMCHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nMCHitsV", nMCHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nMCHitsW", nMCHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nSharedHitsU", 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nSharedHitsV", 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nSharedHitsW", 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nContributingHitsU", 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nContributingHitsV", 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), mcTreeName.c_str(), "nContributingHitsW", 0));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), mcTreeName.c_str()));
        }

        ++mcIndex;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterValidationAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteToTree",
        m_writeToTree));
    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    }

    std::string criteria;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Criteria", criteria));
    std::transform(criteria.begin(), criteria.end(), criteria.begin(), [](unsigned char c){ return std::tolower(c); });
    if (!criteria.empty())
    {
        if (criteria == "neutrino")
        {
            m_criteria = LArMCParticleHelper::IsBeamNeutrinoFinalState;
        }
        else if (criteria == "cosmic")
        {
            m_criteria = LArMCParticleHelper::IsCosmicRay;
        }
        else if (criteria == "testbeam")
        {
            m_criteria = LArMCParticleHelper::IsBeamParticle;
        }
        else
        {
            std::cout << "ClusterValidationAlgorithm: Error. Criteria \'" << criteria << "\'not recognised" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

