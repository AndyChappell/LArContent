/**
 *  @file   larpandoracontent/LArUtility/TpcVolumeIteratorAlgorithm.cc
 *
 *  @brief  Implementation of the TPC volume iterator algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArTpcGeometryHelper.h"
#include "larpandoracontent/LArUtility/TpcVolumeIteratorAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TpcVolumeIteratorAlgorithm::TpcVolumeIteratorAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TpcVolumeIteratorAlgorithm::Run()
{
    // If specified, change the current calo hit list, i.e. the input to the clustering algorithm
    std::string originalCaloHitListName;

    for (const std::string &view : m_viewVector)
    {
        for (const int cryo : m_cryoStatVector)
        {
            for (const int tpc : m_tpcVector)
            {
                for (const int child : m_tpcChildVolumeVector)
                {
                    const std::string suffix{"_" + std::to_string(cryo) + "_" + std::to_string(tpc) + "_" + std::to_string(child)};
                    const std::string caloHitListName{m_caloHitListPrefix + view + suffix };

                    if (STATUS_CODE_NOT_FOUND == PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, caloHitListName))
                        continue;
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_childAlgorithmName));
                    const std::string clusterListName{m_clusterListPrefix + view + suffix };
                    PandoraContentApi::SaveList<Cluster>(*this, m_intermediateClusterListName, clusterListName);
                }
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TpcVolumeIteratorAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "ClusteringWrapper", m_childAlgorithmName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListPrefix", m_caloHitListPrefix));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "IntermediateClusterListName", m_intermediateClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListPrefix", m_clusterListPrefix));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ViewList", m_viewVector));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CryoStatList", m_cryoStatVector));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "TpcList", m_tpcVector));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "TpcChildVolumeList", m_tpcChildVolumeVector));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

