/**
 *  @file   larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

MCMonitoringAlgorithm::MCMonitoringAlgorithm() :
    m_caloHitListName(""),
    m_mcListName("Input"),
    m_visualise(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCMonitoringAlgorithm::~MCMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::Run()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    }

    this->BuildMCHitMap();

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::BuildMCHitMap()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcListName, pMCParticleList));

    if (!pMCParticleList || pMCParticleList->empty())
        return STATUS_CODE_SUCCESS;

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    m_mcHitsMap.clear();
    for (const MCParticle *const pMC : *pMCParticleList)
        m_mcHitsMap[pMC] = CaloHitList();

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            m_mcHitsMap[pMC].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
            continue;
        }
    }

    for (const auto &[pMC, hits] : m_mcHitsMap)
    {
        const MCParticleList &parentList{pMC->GetParentList()};
        std::string desc{""};
        if (parentList.size() == 1)
        {
            std::ostringstream oss;
            oss << parentList.front();
            desc += std::to_string(parentList.front()->GetParticleId()) + " (" + oss.str() + ")" + " -> ";
        }
        else
        {
            if (std::abs(pMC->GetParticleId()) == 14)
            {
                std::cout << "Nu:" << std::endl;
                this->PrintTree(pMC, 0);
                std::cout << "==========================" << std::endl;
            }
        }
        desc += std::to_string(pMC->GetParticleId());

        if (m_visualise && !hits.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, desc, AUTOITER));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCMonitoringAlgorithm::PrintTree(const MCParticle *const pMC, const int depth)
{
    std::string tab(2 * depth, ' ');
    std::cout << tab << pMC->GetParticleId() << ": " << m_mcHitsMap[pMC].size() << std::endl;
    for (const MCParticle *const pChild : pMC->GetDaughterList())
        this->PrintTree(pChild, depth + 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCListName", m_mcListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
