/**
 *  @file   larpandoracontent/LArCheating/CheatingPartialClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingPartialClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingPartialClusterCreationAlgorithm::CheatingPartialClusterCreationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPartialClusterCreationAlgorithm::Run()
{
    MCParticleToHitListMap mcParticleToHitListMap;
    this->GetMCParticleToHitListMap(mcParticleToHitListMap);
    this->CreateClusters(mcParticleToHitListMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPartialClusterCreationAlgorithm::GetMCParticleToHitListMap(MCParticleToHitListMap &mcParticleToHitListMap) const
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            mcParticleToHitListMap[pMCParticle].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPartialClusterCreationAlgorithm::CreateClusters(const MCParticleToHitListMap &mcParticleToHitListMap) const
{
    MCParticleList mcParticleList;
    for (const auto & [ pMCParticle, caloHitList ] : mcParticleToHitListMap)
    {
        // Skip over MC particles with a PDG code in the exclusion list
        if (std::find(m_pdgExclusionVector.begin(), m_pdgExclusionVector.end(), pMCParticle->GetParticleId()) != m_pdgExclusionVector.end())
            continue;
        if (caloHitList.empty())
            continue;

        const Cluster *pCluster{nullptr};
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList = caloHitList;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPartialClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PdgExclusionList", m_pdgExclusionVector));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
