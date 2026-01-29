/**
 *  @file   larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the secondary validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterValidationAlgorithm::ClusterValidationAlgorithm() :
    m_eventNumber(-1),
    m_writeFile(true),
    m_fileName("HierarchyRecoPerformance.root"),
    m_treeName("tree"),
    m_minRecoHits{30},
    m_minRecoHitsPerView{10},
    m_minRecoGoodViews{2},
    m_removeRecoNeutrons{true}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterValidationAlgorithm::~ClusterValidationAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterValidationAlgorithm::Run()
{
    std::string m_caloHitListName2D("CaloHitList2D");
    std::string m_mcParticleListName("Input");

    ++m_eventNumber;

    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName2D, pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    const VertexList *pNuVertexList(nullptr);
    PandoraContentApi::GetList(*this, "NeutrinoVertices3D", pNuVertexList);

    // Get nu vertex accuracy (for later tree filling)
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    if (mcNeutrinoVector.size() != 1)
        return STATUS_CODE_SUCCESS;

    const MCParticle *const pMCNu(mcNeutrinoVector.front());

    float nuVertexAccuracy(-1.f);

    if (pNuVertexList && !pNuVertexList->empty())
    {
        const CartesianVector &trueNuVertex(pMCNu->GetVertex());
        const CartesianVector &recoNuVertex(pNuVertexList->front()->GetPosition());

        nuVertexAccuracy = (trueNuVertex - recoNuVertex).GetMagnitude();
    }

    // Get target MCParticles
    LArHierarchyHelper::FoldingParameters foldParameters;
    foldParameters.m_foldToLeadingShowers = true;

    const LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria recoCriteria(
        m_minRecoHits, m_minRecoHitsPerView, m_minRecoGoodViews, m_removeRecoNeutrons);

    LArHierarchyHelper::MCHierarchy mcHierarchy(recoCriteria);
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList2D, foldParameters, mcHierarchy);

    // Get MC Neutrino - e.g. neutrinos
    MCParticleList rootMCParticles;
    mcHierarchy.GetRootMCParticles(rootMCParticles);

    // Collect all MCNodes
    for (const MCParticle *const pRootParticle : rootMCParticles)
    {
        // Get all MCNodes in this hierarchy
        LArHierarchyHelper::MCHierarchy::NodeVector hierarchyNodes;
        mcHierarchy.GetFlattenedNodes(pRootParticle, hierarchyNodes);

        GenerationMap generationMap;
        this->BuildVisibleHierarchy(pRootParticle, hierarchyNodes, 1, generationMap);

        // /////////////////////////////////
        // for (auto entry : generationMap)
        // {
        //     std::cout << "(TrackId, PDG): (" <<  ((size_t)(intptr_t *)entry.first->GetUid()) << ", " << entry.first->GetParticleId() << "), Tier:" << entry.second << std::endl;
        // }
        // /////////////////////////////////

        ClusterMatchMap clusterMatchMapU, clusterMatchMapV, clusterMatchMapW;
        this->MatchViewClusters(generationMap, hierarchyNodes, TPC_VIEW_U, clusterMatchMapU);
        this->MatchViewClusters(generationMap, hierarchyNodes, TPC_VIEW_V, clusterMatchMapV);
        this->MatchViewClusters(generationMap, hierarchyNodes, TPC_VIEW_W, clusterMatchMapW);

        // /////////////////////////////////
        // std::cout << "U Matches: " << std::endl;
        // for (auto entry : clusterMatchMapU)
        // {
        //     std::cout << "(TrackId, PDG): (" <<  ((size_t)(intptr_t *)entry.first->GetUid()) << ", " << entry.first->GetParticleId() << "), nMatches:" << entry.second.size() << std::endl;
        // }

        // std::cout << "V Matches: " << std::endl;
        // for (auto entry : clusterMatchMapV)
        // {
        //     std::cout << "(TrackId, PDG): (" <<  ((size_t)(intptr_t *)entry.first->GetUid()) << ", " << entry.first->GetParticleId() << "), nMatches:" << entry.second.size() << std::endl;
        // }

        // std::cout << "W Matches: " << std::endl;
        // for (auto entry : clusterMatchMapW)
        // {
        //     std::cout << "(TrackId, PDG): (" <<  ((size_t)(intptr_t *)entry.first->GetUid()) << ", " << entry.first->GetParticleId() << "), nMatches:" << entry.second.size() << std::endl;
        // }
        // /////////////////////////////////

        this->FillTree(generationMap, hierarchyNodes, clusterMatchMapU, clusterMatchMapV, clusterMatchMapW, nuVertexAccuracy, pMCNu);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::BuildVisibleHierarchy(const MCParticle *const pMCParent, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes,
    const int childTier, GenerationMap &generationMap)
{
    for (const MCParticle *const pMCChild : pMCParent->GetDaughterList())
    {        
        // Find MC node
        const LArHierarchyHelper::MCHierarchy::Node *pChildMCNode(nullptr);
        this->FindMCNode(pMCChild, hierarchyNodes, pChildMCNode);
        
        // If we can't find a node or if it isn't reco'able skip it
        if ((pChildMCNode == nullptr) || (!pChildMCNode->IsReconstructable()))
        {   
            this->BuildVisibleHierarchy(pMCChild, hierarchyNodes, childTier, generationMap);
        }
        else
        {
            generationMap[pMCChild] = childTier;
            this->BuildVisibleHierarchy(pMCChild, hierarchyNodes, (childTier + 1), generationMap);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::MatchViewClusters(const GenerationMap &generationMap, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes, 
    const HitType hitType, ClusterMatchMap &clusterMatchMap)
{
    std::string m_clusterListNameU("ClustersU");
    std::string m_clusterListNameV("ClustersV");
    std::string m_clusterListNameW("ClustersW");
    const std::string clusterListName(hitType == TPC_VIEW_U ? m_clusterListNameU : hitType == TPC_VIEW_V ? m_clusterListNameV : m_clusterListNameW);

    for (auto &entry : generationMap)
        clusterMatchMap[entry.first] = std::vector<std::pair<const pandora::Cluster*, int>>();

    const ClusterList *pClusterList(nullptr);
    if (PandoraContentApi::GetList(*this, clusterListName, pClusterList) != STATUS_CODE_SUCCESS)
        return;

    for (const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList clusterHitList;
        LArClusterHelper::GetAllHits(pCluster, clusterHitList);

        int bestSharedHits(0); float bestSharedEnergy(-1.f); const MCParticle *pBestMatch(nullptr);

        for (const LArHierarchyHelper::MCHierarchy::Node *const pMCNode : hierarchyNodes)
        {
            CaloHitList targetHitList;
            const MCParticle *const pNodeMCParticle(pMCNode->GetMCParticles().front());
            this->GetHitsOfType(pMCNode->GetCaloHits(), hitType, targetHitList);;

            // Tried to user the set_intersection function but it didnt like the ordering?
            CaloHitVector intersection;
            for (const CaloHit *const pCaloHit : clusterHitList)
            {
                if (std::find(targetHitList.begin(), targetHitList.end(), pCaloHit) != targetHitList.end())
                    intersection.push_back(pCaloHit);
            }

            const int nSharedHits(intersection.size());

            if (nSharedHits == 0)
            {
                continue;
            }
            else if (nSharedHits == bestSharedHits)
            {
                float totalEnergy(0.f);
                for (const CaloHit *const pCaloHit : intersection)
                    totalEnergy += pCaloHit->GetElectromagneticEnergy();

                if (totalEnergy > bestSharedEnergy)
                {
                    bestSharedEnergy = totalEnergy;
                    pBestMatch = pNodeMCParticle;
                }
            }
            else if (nSharedHits > bestSharedHits)
            {
                float totalEnergy(0.f);
                for (const CaloHit *const pCaloHit : intersection)
                    totalEnergy += pCaloHit->GetElectromagneticEnergy();

                bestSharedHits = nSharedHits;
                bestSharedEnergy = totalEnergy;
                pBestMatch = pNodeMCParticle;
            }
        }

        if (pBestMatch)
            clusterMatchMap[pBestMatch].push_back(std::make_pair(pCluster, bestSharedHits));
    }

    // sort vector so that best match is first!
    for (auto &entry : clusterMatchMap)
    {
        std::sort(entry.second.begin(), entry.second.end(),
                  [](std::pair<const Cluster*, int> lhs, std::pair<const Cluster*, int> rhs)
                  { return lhs.second > rhs.second; });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::FillTree(const GenerationMap &generationMap, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes, 
    const ClusterMatchMap &clusterMatchMapU, const ClusterMatchMap &clusterMatchMapV, const ClusterMatchMap &clusterMatchMapW, 
    const float nuVertexAccuracy, const MCParticle *const pMCNu)
{
    int nTargets(generationMap.size());

    for (auto &entry : generationMap)
    {
        IntVector nSharedHitsU, nSharedHitsV, nSharedHitsW;
        FloatVector completenessU, completenessV, completenessW;
        FloatVector purityU, purityV, purityW;
        IntVector isViewReconstructable;

        const MCParticle *const pMCTarget(entry.first);

        const LArHierarchyHelper::MCHierarchy::Node *pMCNode(nullptr);
        this->FindMCNode(entry.first, hierarchyNodes, pMCNode);

        CaloHitList targetHitList(pMCNode->GetCaloHits());
        CaloHitList targetHitListU, targetHitListV, targetHitListW;
        this->GetHitsOfType(targetHitList, TPC_VIEW_U, targetHitListU);
        this->GetHitsOfType(targetHitList, TPC_VIEW_V, targetHitListV);
        this->GetHitsOfType(targetHitList, TPC_VIEW_W, targetHitListW);

        std::vector<std::pair<const Cluster*, int>> matchesU(clusterMatchMapU.at(pMCTarget));
        std::vector<std::pair<const Cluster*, int>> matchesV(clusterMatchMapV.at(pMCTarget));
        std::vector<std::pair<const Cluster*, int>> matchesW(clusterMatchMapW.at(pMCTarget));

        this->GetMatchingMetrics(targetHitListU, matchesU, completenessU, purityU, nSharedHitsU, isViewReconstructable);
        this->GetMatchingMetrics(targetHitListV, matchesV, completenessV, purityV, nSharedHitsV, isViewReconstructable);
        this->GetMatchingMetrics(targetHitListW, matchesW, completenessW, purityW, nSharedHitsW, isViewReconstructable);

        int nMCHitsU(targetHitListU.size()), nMCHitsV(targetHitListV.size()), nMCHitsW(targetHitListW.size());

        // /////////////////////////////////
        // std::cout << "-------------------" << std::endl;
        //     std::cout << "(TrackId, PDG): (" <<  ((size_t)(intptr_t *)pMCTarget->GetUid()) << ", " << pMCTarget->GetParticleId() << "), Tier:" << entry.second << std::endl;
        // std::cout << "-------------------" << std::endl;
        // for (unsigned int i = 0; i < matchesU.size(); ++i)
        // {
        //     std::cout << "(nShared, completeness, purity) = (" << nSharedHitsU.at(i) << ", " << completenessU.at(i) << ", " << purityU.at(i) << ")" << std::endl;
        // }
        // /////////////////////////////////

        // Fill tree
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_eventNumber));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NuPDG", pMCNu->GetParticleId()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NuEnergy", pMCNu->GetEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NTargets", nTargets));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NuVertexAccuracy", nuVertexAccuracy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PDG", pMCTarget->GetParticleId()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Energy", pMCTarget->GetEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NMCHitsU", nMCHitsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NCMHitsV", nMCHitsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NMCHitsW", nMCHitsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VisibleTier", entry.second));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsViewReconstructable", &isViewReconstructable));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NSharedHitsU", &nSharedHitsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NSharedHitsV", &nSharedHitsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NSharedHitsW", &nSharedHitsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "CompletenessU", &completenessU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "CompletenessV", &completenessV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "CompletenessW", &completenessW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PurityU", &purityU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PurityV", &purityV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PurityW", &purityW));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::GetHitsOfType(const CaloHitList &inputList, const HitType hitType, CaloHitList &outputList)
{
    for (const CaloHit *const pCaloHit : inputList)
    {
        if (pCaloHit->GetHitType() == hitType)
            outputList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterValidationAlgorithm::GetMatchingMetrics(const CaloHitList &targetHitList, const std::vector<std::pair<const Cluster*, int>> matches,
    FloatVector &completeness, FloatVector &purity, IntVector &nSharedHits, IntVector &isReconstructableView)
{
    isReconstructableView.push_back((targetHitList.size() < m_minRecoHitsPerView) ? 0 : 1);

    for (auto &match : matches)
    {
        CaloHitList clusterHitList;
        LArClusterHelper::GetAllHits(match.first, clusterHitList);

        nSharedHits.push_back(match.second);
        completeness.push_back(targetHitList.empty() ? 0 : static_cast<float>(match.second) / static_cast<float>(targetHitList.size()));
        purity.push_back(clusterHitList.empty() ? 0 : static_cast<float>(match.second) / static_cast<float>(clusterHitList.size()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Find functions
//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterValidationAlgorithm::FindMCNode(const MCParticle *const pMCParticle, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyMCNodes, 
    const LArHierarchyHelper::MCHierarchy::Node *&pMCNode)
{
    for (const LArHierarchyHelper::MCHierarchy::Node *const pNode : hierarchyMCNodes)
    {
        if (pNode->GetMCParticles().front() == pMCParticle)
        {
            pMCNode = pNode;
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHits", m_minRecoHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHitsPerView", m_minRecoHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoGoodViews", m_minRecoGoodViews));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RemoveRecoNeutrons", m_removeRecoNeutrons));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
