/**
 *  @file   larpandoracontent/LArCheating/CheatingVertexHitAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArEventTopology.h"

#include "larpandoracontent/LArCheating/CheatingVertexHitAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingVertexHitAssociationAlgorithm::CheatingVertexHitAssociationAlgorithm() :
    m_inputCaloHitListName("CaloHitList2D")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexHitAssociationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            mcToHitsMap[pMCParticle].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }

    const ClusterList *pOutputClusterList{nullptr};
    std::string originalClusterListName, tempListName{"temp"};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputClusterList, tempListName));
    for (auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const Cluster *pCluster{nullptr};
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList.emplace_back(caloHitList.front());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

        caloHitList.pop_front();
        for (const CaloHit *const pAssociatedCaloHit : caloHitList)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pAssociatedCaloHit));
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, originalClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalClusterListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingVertexHitAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_inputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
