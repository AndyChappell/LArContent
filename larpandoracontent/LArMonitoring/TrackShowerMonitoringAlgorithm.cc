/**
 *  @file   larpandoracontent/LArMonitoring/TrackShowerMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the track shower monitoring algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TrackShowerMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

//-----------------------------------------------------------------------------------------------------------------------------------------

TrackShowerMonitoringAlgorithm::TrackShowerMonitoringAlgorithm() :
    m_showTruth{false},
    m_showNetworkClass{false}
{
    this->m_viewToNameMap.insert(std::make_pair(TPC_VIEW_U, "U"));
    this->m_viewToNameMap.insert(std::make_pair(TPC_VIEW_V, "V"));
    this->m_viewToNameMap.insert(std::make_pair(TPC_VIEW_W, "W"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true,
                DETECTOR_VIEW_XZ,  -1.f, 1.f, 1.f));

    // Show calo hit truth
    if (m_showTruth)
    {
        this->VisualizeCaloHitTruth();
    }

    // Show specified lists of pfo
    for (const std::string &listName : m_pfoListNames)
    {
        this->VisualizePfoList(listName);
    }

    // Show specified lists of clusters
    for (const std::string &listName : m_clusterListNames)
    {
        this->VisualizeAvailableClusterList(listName);
    }

    // Show calo hit network classification
    if (m_showTruth)
    {
        this->VisualizeNetworkClassification();
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::VisualizePfoList(const std::string &listName) const
{
    const PfoList *pPfoList = nullptr;
    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pPfoList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "TrackShowerMonitoringAlgorithm: current pfo list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pPfoList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "TrackShowerMonitoringAlgorithm: pfo list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    std::map<HitType, ClusterList*> trackClusterLists;
    std::map<HitType, ClusterList*> showerClusterLists;
    const auto hitEnums = {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};
    for (const auto hitType : hitEnums)
    {
        trackClusterLists.insert(std::make_pair(hitType, new ClusterList()));
        showerClusterLists.insert(std::make_pair(hitType, new ClusterList()));
    }

    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        const ClusterList &pfoClusterList = pPfo->GetClusterList();
        if (LArPfoHelper::IsTrack(pPfo))
        {
            for (const Cluster *pCluster : pfoClusterList)
                trackClusterLists.at(LArClusterHelper::GetClusterHitType(pCluster))->push_back(pCluster);
        }
        else if (LArPfoHelper::IsShower(pPfo))
        {
            for (const Cluster *pCluster : pfoClusterList)
                showerClusterLists.at(LArClusterHelper::GetClusterHitType(pCluster))->push_back(pCluster);
        }
    }

    for (const auto [ hitType, clusterList ] : trackClusterLists)
    {
        if (clusterList->size() > 0)
        {
            std::string name = "TrackClusters" + this->m_viewToNameMap.at(hitType);
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), clusterList, name, CYAN));
        }
    }
    for (const auto [ hitType, clusterList ] : showerClusterLists)
    {
        if (clusterList->size() > 0)
        {
            std::string name = "ShowerClusters" + this->m_viewToNameMap.at(hitType);
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), clusterList, name, ORANGE));
        }
    }

    for (const auto hitType : hitEnums)
    {
        delete trackClusterLists.at(hitType);
        delete showerClusterLists.at(hitType);
    }
    trackClusterLists.clear();
    showerClusterLists.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::VisualizeAvailableClusterList(const std::string &listName) const
{
    const ClusterList *pClusterList = nullptr;
    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pClusterList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "TrackShowerMonitoringAlgorithm: current cluster list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pClusterList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "TrackShowerMonitoringAlgorithm: cluster list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    ClusterList availableClusterList;

    for (const Cluster *pCluster : *pClusterList)
    {
        if (pCluster->IsAvailable())
            availableClusterList.push_back(pCluster);
    }

    if (availableClusterList.size() > 0)
    {
        std::string name = "AvailableClusters" + this->m_viewToNameMap.at(LArClusterHelper::GetClusterHitType(
                    availableClusterList.front()));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &availableClusterList, name, GREEN));
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::VisualizeCaloHitTruth() const
{
    const CaloHitList *pCaloHitList2D(nullptr);
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, this->m_caloHitList2DName, pCaloHitList2D))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TrackShowerMonitoringAlgorithm: calo hit list " << this->m_caloHitList2DName << " unavailable." << std::endl;
        return;
    }
    const MCParticleList *pMCParticleList(nullptr);
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pMCParticleList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TrackShowerMonitoringAlgorithm: current mc particle list unavailable." << std::endl;
        return;
    }

    LArMCParticleHelper::PrimaryParameters parameters;
    // Turn off max photo propagation for now, only care about killing off daughters of neutrons
    parameters.m_maxPhotonPropagation = 100.;
    parameters.m_foldBackHierarchy = false;
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
            LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    for (const std::string listName : this->m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList = nullptr;
        PandoraContentApi::GetList(*this, listName, pCaloHitList);
        CaloHitList trackHitsList, showerHitsList, michelHitsList;
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                // Throw away non-reconstructable hits
                if (targetMCParticleToHitsMap.find(pMCParticle) == targetMCParticleToHitsMap.end())
                    continue;
                
                const int pdg = pMCParticle->GetParticleId();
                if (LArMCParticleHelper::IsDescendentOf(pMCParticle, 2112))
                {
                    continue;
                }
                else if (LArMCParticleHelper::IsDescendentOf(pMCParticle, 111))
                {
                    showerHitsList.push_back(pCaloHit);
                }
                else if (std::abs(pdg) == 11)
                {
                    if (LArMCParticleHelper::IsDescendentOf(pMCParticle, 13))
                        michelHitsList.push_back(pCaloHit);
                    else
                        showerHitsList.push_back(pCaloHit);
                }
                else if (std::abs(pdg) == 22)
                {
                    showerHitsList.push_back(pCaloHit);
                }
                else if (std::abs(pdg) == 2212 || std::abs(pdg) > 1e9)
                {
                    trackHitsList.push_back(pCaloHit);
                }
                else if (std::abs(pdg) == 13 || std::abs(pdg) == 211)
                {
                    trackHitsList.push_back(pCaloHit);
                }
                else
                {
                    trackHitsList.push_back(pCaloHit);
                }
            }
            catch (...)
            {
                continue;
            }
        }
        if (trackHitsList.size() > 0)
        {
            std::string tracksName = "MCTracks" + this->m_viewToNameMap.at(trackHitsList.front()->GetHitType());
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHitsList, tracksName, BLUE));
        }
        if (showerHitsList.size() > 0)
        {
            std::string showersName = "MCShowers" + this->m_viewToNameMap.at(showerHitsList.front()->GetHitType());
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHitsList, showersName, RED));
        }
        if (michelHitsList.size() > 0)
        {
            std::string michelsName = "MCMichels" + this->m_viewToNameMap.at(michelHitsList.front()->GetHitType());
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &michelHitsList, michelsName, YELLOW));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::VisualizeNetworkClassification() const
{
    for (const std::string listName : this->m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList = nullptr;
        PandoraContentApi::GetList(*this, listName, pCaloHitList);
        CaloHitList trackHitsList, showerHitsList;
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            auto properties = pCaloHit->GetPropertiesMap();
            if (properties["Ptrack"] >= properties["Pshower"])
                trackHitsList.push_back(pCaloHit);
            else
                showerHitsList.push_back(pCaloHit);
        }
        if (trackHitsList.size() > 0)
        {
            std::string tracksName = "NetworkTracks" + this->m_viewToNameMap.at(trackHitsList.front()->GetHitType());
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHitsList, tracksName, MAGENTA));
        }
        if (showerHitsList.size() > 0)
        {
            std::string showersName = "NetworkShowers" + this->m_viewToNameMap.at(showerHitsList.front()->GetHitType());
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHitsList, showersName, PINK));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowTruth", m_showTruth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowNetworkClass", m_showNetworkClass));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitList2DName", m_caloHitList2DName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "PfoListNames", m_pfoListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content