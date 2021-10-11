/**
 *  @file   larpandoracontent/LArMonitoring/AtmosphericsMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/AtmosphericsMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

AtmosphericsMonitoringAlgorithm::AtmosphericsMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AtmosphericsMonitoringAlgorithm::~AtmosphericsMonitoringAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_outputTreeName, m_outputFileName, "RECREATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "AtmosphericMonitoringAlgorithm: Unable to write to ROOT tree" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AtmosphericsMonitoringAlgorithm::Run()
{
    this->ConstructRootTree();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AtmosphericsMonitoringAlgorithm::ConstructRootTree() const
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
    if (!primaries.empty())
    {
        const MCParticle *primary{primaries.front()};
        const MCParticleList &parents{primary->GetParentList()};
        if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
            pTrueNeutrino = parents.front();
    }

    if (pTrueNeutrino)
    {
        const int flavour{pTrueNeutrino->GetParticleId()};
        const float energy{pTrueNeutrino->GetEnergy()};
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "flavour", flavour));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "energy", energy));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_outputTreeName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AtmosphericsMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFileName", m_outputFileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTreeName", m_outputTreeName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
