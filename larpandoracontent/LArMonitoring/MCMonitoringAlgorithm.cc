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
    m_visualise(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCMonitoringAlgorithm::~MCMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::Run()
{

    this->BuildMCHitMap();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCMonitoringAlgorithm::BuildMCHitMap()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    }

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

    MCHitsMap pizeroHitsMap, otherHitsMap;
    for (const auto &[pMC, hits] : m_mcHitsMap)
    {
        if (pMC->GetParticleId() == 111)
            continue;
        MCParticleList parentList{pMC->GetParentList()};
        bool found_pizero{false};
        while (!parentList.empty())
        {
            if (parentList.front()->GetParticleId() == 111)
            {
                found_pizero = true;
                if (pizeroHitsMap.find(parentList.front()) != pizeroHitsMap.end())
                {
                    pizeroHitsMap[parentList.front()].insert(pizeroHitsMap[parentList.front()].end(),
                        hits.begin(), hits.end());
                }
                else
                {
                    pizeroHitsMap[parentList.front()] = CaloHitList();
                    pizeroHitsMap[parentList.front()].insert(pizeroHitsMap[parentList.front()].end(),
                        hits.begin(), hits.end());
                }
                break;
            }
            parentList = parentList.front()->GetParentList();
        }
        if (!found_pizero)
        {
            otherHitsMap[pMC] = hits;
        }
    }
    for (const auto &[pMC, hits] : pizeroHitsMap)
    {
        if (hits.empty())
            continue;

        if (m_visualise)
        {
            std::ostringstream oss;
            oss << pMC->GetParticleId() << " (" << hits.size() << ")";
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, oss.str(), AUTOITER));
        }

        std::cout << "Found visible pi-zero with " << hits.size() << " hits and momentum "
                  << pMC->GetMomentum().GetMagnitude() << " GeV." << std::endl;
    }
    for (const auto &[pMC, hits] : otherHitsMap)
    {
        if (hits.empty())
            continue;

        if (m_visualise)
        {
            std::ostringstream oss;
            oss << pMC->GetParticleId() << " (" << hits.size() << ")";
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, oss.str(), BLACK));
        }
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
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
            if (parentList.front()->GetParticleId() == 111)
            {
                std::cout << "MCParticle: PDG " << pMC->GetParticleId() << " with parent PDG 111 has " << hits.size() << " hits." << std::endl;
            }
        }
        desc += std::to_string(pMC->GetParticleId());

        if (m_visualise && !hits.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, desc, AUTOITER));
        }
    }

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
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
