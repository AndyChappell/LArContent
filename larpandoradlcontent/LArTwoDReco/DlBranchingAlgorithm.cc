/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlBranchingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArTwoDReco/DlBranchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlBranchingAlgorithm::DlBranchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlBranchingAlgorithm::~DlBranchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlBranchingAlgorithm::Run()
{
    const ClusterList *pTrackClusterList{nullptr};
    // Set the input track cluster list as current
    PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_trackClusterListName);
    StatusCode code{PandoraContentApi::GetCurrentList(*this, pTrackClusterList)};
    if (code == STATUS_CODE_SUCCESS)
    {
        for (const auto &alg : m_trackAlgorithms)
        {   // ATTN - The clustering algorithms replace the current cluster list as they go
            PandoraContentApi::GetCurrentList(*this, pTrackClusterList);
            if (!pTrackClusterList->empty())
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, alg));
            }
        }
        // Save the current cluster list to the target output list
        PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName);
    }
    else if (code != STATUS_CODE_NOT_INITIALIZED)
    {
        return code;
    }

    const ClusterList *pShowerClusterList{nullptr};
    // Set the input shower cluster list as current
    PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_showerClusterListName);
    code = PandoraContentApi::GetCurrentList(*this, pShowerClusterList);
    if (code == STATUS_CODE_SUCCESS)
    {
        for (const auto &alg : m_showerAlgorithms)
        {
            PandoraContentApi::GetCurrentList(*this, pShowerClusterList);
            if (!pShowerClusterList->empty())
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, alg));
            }
        }
        // Append the current cluster list to the target output list
        PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName);
    }
    else if (code != STATUS_CODE_NOT_INITIALIZED)
    {
        return code;
    }

    // Set the current list to the final output list
    PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputClusterListName);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlBranchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackClusterListName", m_trackClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerClusterListName", m_showerClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "TrackAlgorithms", m_trackAlgorithms));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "ShowerAlgorithms", m_showerAlgorithms));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
