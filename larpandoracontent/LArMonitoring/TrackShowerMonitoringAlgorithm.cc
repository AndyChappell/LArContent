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

//-----------------------------------------------------------------------------------------------------------------------------------------

TrackShowerMonitoringAlgorithm::~TrackShowerMonitoringAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "net_tree", "net.root", "UPDATE"));
    }
    catch(const StatusCodeException&)
    {
        std::cout << "TrackShowerMonitoringAlgorithm: Unable to write tree net_tree to file net.root" << std::endl;
    }
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
        // Create classification histograms
        this->SerializePfoClassification(listName);
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
    try
    {
        this->GetPfoList(listName, pPfoList);
        std::cout << "Found pfo list " << pPfoList << std::endl;
    }
    catch(const StatusCodeException&)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TrackShowerMonitoringAlgorithm: pfo list \'" << listName << "\' unavailable." << std::endl;
        return;
    }

    ViewToClusterListMap trackClusterLists, showerClusterLists;
    MakeClusterMaps(trackClusterLists, showerClusterLists);

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

    DestroyClusterMaps(trackClusterLists, showerClusterLists);
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
                Classification cls{this->GetTruthTag(*pCaloHit, targetMCParticleToHitsMap)};
                // Throw away non-reconstructable hits
                if (cls == NON_RECO)
                    continue;
                else if (cls == TRACK)
                    trackHitsList.push_back(pCaloHit);
                else if (cls == SHOWER)
                    showerHitsList.push_back(pCaloHit);
                else if (cls == MICHEL)
                    michelHitsList.push_back(pCaloHit);
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
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHitsList, tracksName, DARKBLUE));
        }
        if (showerHitsList.size() > 0)
        {
            std::string showersName = "NetworkShowers" + this->m_viewToNameMap.at(showerHitsList.front()->GetHitType());
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHitsList, showersName, DARKRED));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::SerializePfoClassification(const std::string &listName) const
{
    static int e{0};
    const PfoList *pPfoList = nullptr;
    try
    {
        this->GetPfoList(listName, pPfoList);
    }
    catch(const StatusCodeException&)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TrackShowerMonitoringAlgorithm: pfo list \'" << listName << "\' unavailable." << std::endl;
        e++;
        return;
    }

    const CaloHitList *pCaloHitList2D(nullptr);
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, this->m_caloHitList2DName, pCaloHitList2D))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TrackShowerMonitoringAlgorithm: calo hit list " << this->m_caloHitList2DName << " unavailable." << std::endl;
        e++;
        return;
    }
    const MCParticleList *pMCParticleList(nullptr);
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pMCParticleList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "TrackShowerMonitoringAlgorithm: current mc particle list unavailable." << std::endl;
        e++;
        return;
    }

    LArMCParticleHelper::PrimaryParameters parameters;
    // Turn off max photo propagation for now, only care about killing off daughters of neutrons
    parameters.m_maxPhotonPropagation = 100.;
    parameters.m_foldBackHierarchy = false;
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
            LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    int p{0};
    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        const ClusterList &pfoClusterList = pPfo->GetClusterList();
        if (LArPfoHelper::IsTrack(pPfo))
        {
            int c{0}, pfoTrack{1};
            for (const Cluster *pCluster : pfoClusterList)
            {
                FloatVector clusterProb, pfoProb;
                //if (LArClusterHelper::GetClusterHitType(pCluster) == TPC_VIEW_W)
                {
                    float trackEnergy{0.f}, clusterEnergy{0.f};
                    CaloHitList caloHitList;
                    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);;
                    for (const CaloHit *pCaloHit : caloHitList)
                    {
                        Classification cls{this->GetTruthTag(*pCaloHit, targetMCParticleToHitsMap)};
                        float energy = pCaloHit->GetInputEnergy();
                        energy = energy > 0.f ? energy : 0.f;
                        if (cls == TRACK)
                        {
                            trackEnergy += energy;
                            clusterEnergy += energy;
                        }
                        else
                        {
                            clusterEnergy += energy;
                        }
                        auto properties = pCaloHit->GetPropertiesMap();
                        clusterProb.push_back(properties["Ptrack"]);
                        pfoProb.push_back(properties["Ptrack"]);
                    }
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "net_tree", "eventId", e));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "net_tree", "pfoId", p));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "net_tree", "clusterId", c));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "net_tree", "pfoIsTrack", pfoTrack));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "net_tree", "clusterTrackEnergy", trackEnergy));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "net_tree", "clusterEnergy", clusterEnergy));
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "net_tree", "clusterProb", &clusterProb));
                    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "net_tree"));
                    c++;
                }
            }
        }
        else if (LArPfoHelper::IsShower(pPfo))
        {
            // TODO
        }
        p++;
    }
    e++;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::GetPfoList(const std::string &listName, const PfoList* &pPfoList) const
{
    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pPfoList))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pPfoList))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

Classification TrackShowerMonitoringAlgorithm::GetTruthTag(const pandora::CaloHit &caloHit,
        const LArMCParticleHelper::MCContributionMap &targetMCParticleToHitsMap) const
{
    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(&caloHit));
    // Throw away non-reconstructable hits
    if (targetMCParticleToHitsMap.find(pMCParticle) == targetMCParticleToHitsMap.end())
        return NON_RECO;
    
    const int pdg = pMCParticle->GetParticleId();
    if (LArMCParticleHelper::IsDescendentOf(pMCParticle, 2112))
    {
        return NON_RECO;
    }
    else if (LArMCParticleHelper::IsDescendentOf(pMCParticle, 111))
    {
        return SHOWER;
    }
    else if (std::abs(pdg) == 11)
    {
        if (LArMCParticleHelper::IsDescendentOf(pMCParticle, 13))
            return MICHEL;
        else
            return SHOWER;
    }
    else if (std::abs(pdg) == 22)
    {
        return SHOWER;
    }
    else if (std::abs(pdg) == 2212 || std::abs(pdg) > 1e9)
    {
        return TRACK;
    }
    else if (std::abs(pdg) == 13 || std::abs(pdg) == 211)
    {
        return TRACK;
    }
    else
    {
        return TRACK;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::MakeClusterMaps(ViewToClusterListMap& trackClusterLists, ViewToClusterListMap& showerClusterLists)
    const
{
    const auto hitEnums = {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};
    for (const auto hitType : hitEnums)
    {
        trackClusterLists.insert(std::make_pair(hitType, new ClusterList()));
        showerClusterLists.insert(std::make_pair(hitType, new ClusterList()));
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerMonitoringAlgorithm::DestroyClusterMaps(ViewToClusterListMap& trackClusterLists, ViewToClusterListMap& showerClusterLists)
    const
{
    const auto hitEnums = {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};
    for (const auto hitType : hitEnums)
    {
        delete trackClusterLists.at(hitType);
        delete showerClusterLists.at(hitType);
    }
    trackClusterLists.clear();
    showerClusterLists.clear();
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
