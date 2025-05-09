/**
 *  @file   larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the mc particle monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArTrackPfo.h"

using namespace pandora;

namespace lar_content
{

MCMonitoringAlgorithm::MCMonitoringAlgorithm() :
    m_caloHitListName{"CaloHitList2D"},
    m_mcParticleListName{"Input"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::Run()
{
    std::cout << "Event" << std::endl;
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    for (const MCParticle *const pMC : *pMCParticleList)
        mcToHitsMap[pMC] = CaloHitList();

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const auto mcWeights{pCaloHit->GetMCParticleWeightMap()};
            for (const auto &[pMC, weight] : mcWeights)
                mcToHitsMap.at(pMC).emplace_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
        }
    }

    std::map<const MCParticle*, const MCParticle*> mcToRootMap;
    for (const MCParticle *const pMC : *pMCParticleList)
        mcToRootMap[pMC] = LArMCParticleHelper::GetParentMCParticle(pMC);

    std::set<const MCParticle *> rootNus;
    for (const auto &[pChild, pRoot] : mcToRootMap)
        rootNus.insert(pRoot);
    for (const MCParticle *const pRoot : rootNus)
    {
        const int pdg{std::abs(pRoot->GetParticleId())};
        if (pdg == 12 || pdg == 14 || pdg == 16)
            this->PrintHierarchy(pRoot, "");
    }

    std::cout << "------------------------------------------------------------------------------------------------" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCMonitoringAlgorithm::PrintHierarchy(const MCParticle *const pMC, const std::string &tab) const
{
    std::cout << tab << pMC->GetParticleId() << " " << pMC->GetDaughterList().size() << std::endl;
    for (const MCParticle *const pChild : pMC->GetDaughterList())
        this->PrintHierarchy(pChild, tab + "  ");
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
