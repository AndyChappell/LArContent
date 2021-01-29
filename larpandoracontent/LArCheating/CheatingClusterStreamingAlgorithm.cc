/**
 *  @file   larpandoradlcontent/LArCheating/CheatingClusterStreamingAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster streaming algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingClusterStreamingAlgorithm.h"

#include "Helpers/MCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterStreamingAlgorithm::CheatingClusterStreamingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CheatingClusterStreamingAlgorithm::~CheatingClusterStreamingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterStreamingAlgorithm::Run()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    if (!pMCParticleList || pMCParticleList->empty())
        return STATUS_CODE_NOT_FOUND;

    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    std::string originalClusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalClusterListName));

    ClusterList trackClusterList, showerClusterList;
    for (const Cluster *pCluster : *pClusterList)
    {
        const MCParticle *pMC{nullptr};
        try
        {
            // ATTN: Not all hits have populated MC weight lists, so small clusters may not have matching MC
            pMC = MCParticleHelper::GetMainMCParticle(pCluster);
        }
        catch (const StatusCodeException&)
        {
        }
    
        bool isTrueTrack{true};
        if (pMC)
        {   // Assume clusters without a matching MC are track like
            int pdg{std::abs(pMC->GetParticleId())};
            if (pdg == PHOTON || pdg == E_MINUS)
                isTrueTrack = false;
        }
        if (isTrueTrack)
            trackClusterList.emplace_back(pCluster);
        else
            showerClusterList.emplace_back(pCluster);
    }

    // ATTN - We're ok with saving empty lists here and allowing future algorithms to simply do nothing if there are no clusters
    // Moves the subset of clusters in the cluster list from the old list to the new list
    PandoraContentApi::SaveList<ClusterList>(*this, originalClusterListName, m_trackClusterListName, trackClusterList);
    PandoraContentApi::SaveList<ClusterList>(*this, originalClusterListName, m_showerClusterListName, showerClusterList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterStreamingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackClusterListName", m_trackClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerClusterListName", m_showerClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
