/**
 *  @file   larpandoracontent/LArMonitoring/NeutrinoMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/NeutrinoMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

NeutrinoMonitoringAlgorithm::NeutrinoMonitoringAlgorithm() :
    m_caloHitListName("CaloHitList2D")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoMonitoringAlgorithm::~NeutrinoMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoMonitoringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (pCaloHitList && !pCaloHitList->empty())
    {
        LArMCParticleHelper::MCContributionMap mcHitsMap;
        this->MakeSelection(pCaloHitList, mcHitsMap);
        std::cout << "NeutrinoMonitoringAlgorithm::Run - Found " << mcHitsMap.size() << " neutrino induced particles" << std::endl;
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

        for (const auto &[pMC, caloHits] : mcHitsMap)
        {
            const int pdg{std::abs(pMC->GetParticleId())};
            CaloHitList uHits, vHits, wHits;
            for (const CaloHit *const pCaloHit : caloHits)
            {
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                        uHits.emplace_back(pCaloHit);
                        break;
                    case TPC_VIEW_V:
                        vHits.emplace_back(pCaloHit);
                        break;
                    case TPC_VIEW_W:
                        wHits.emplace_back(pCaloHit);
                        break;
                    default:
                        break;
                }
            }
            if (!uHits.empty())
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits, "u_" + std::to_string(pdg), AUTOITER));
            }
            if (!vHits.empty())
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits, "v_" + std::to_string(pdg), AUTOITER));
            }
            if (!wHits.empty())
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "w_" + std::to_string(pdg), AUTOITER));
            }
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoMonitoringAlgorithm::MakeSelection(const CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcMap) const
{
    int matchedHits(0), unmatchedHits(0);
    CaloHitList matched, unmatched;
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            const MCParticleList &parentList{pMC->GetParentList()};
            if (!parentList.empty())
            {
                const MCParticle *const pParent(parentList.front());
                const int pdg{std::abs(pParent->GetParticleId())};
                if (pdg == 12 || pdg == 14 || pdg == 16)
                    std::cout << "Found a neutrino with pdg " << pdg << std::endl;
            }
            const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMC));
            if (LArMCParticleHelper::IsNeutrino(pParentMCParticle) || LArMCParticleHelper::IsBeamParticle(pParentMCParticle))
            {
                mcMap[pMC].emplace_back(pCaloHit);
                matched.emplace_back(pCaloHit);
            }
            ++matchedHits;
        }
        catch (const StatusCodeException &)
        {
            if (pCaloHit->GetHitType() == TPC_VIEW_W)
                unmatched.emplace_back(pCaloHit);
            ++unmatchedHits;
        }
    }
    std::cout << "NeutrinoMonitoringAlgorithm::MakeSelection - Matched " << matchedHits << " hits, unmatched " << unmatchedHits << std::endl;
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &matched, "matched", RED));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &unmatched, "unmatched", BLACK));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
