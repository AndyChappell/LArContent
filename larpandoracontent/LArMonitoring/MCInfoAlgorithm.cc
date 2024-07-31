/**
 *  @file   larpandoracontent/LArMonitoring/MCInfoAlgorithm.cc
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

#include "larpandoracontent/LArMonitoring/MCInfoAlgorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArTrackPfo.h"

using namespace pandora;

namespace lar_content
{

MCInfoAlgorithm::MCInfoAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCInfoAlgorithm::Run()
{
    std::cout << "---MC-INFO-----------------------------------------------------------------------" << std::endl;
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    const MCParticle *pTrueNeutrino{nullptr};
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            while (!pMC->GetParentList().empty())
            {
                // Check if primary
                if (LArMCParticleHelper::IsNeutrino(pMC->GetParentList().front()))
                    mcToHitsMap[pMC].emplace_back(pCaloHit);
                pMC = pMC->GetParentList().front();
            }

            if (!pTrueNeutrino && LArMCParticleHelper::IsNeutrino(pMC))
                pTrueNeutrino = pMC;
        }
        catch (StatusCodeException &)
        {
        }
    }

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    std::cout << "Nu PDG: " << pTrueNeutrino->GetParticleId() << " E: " << pTrueNeutrino->GetEnergy() << " P:" <<
        pTrueNeutrino->GetMomentum().GetMagnitude() << std::endl;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        int primaryHits{0}, secondaryHits{0};
        CaloHitList uHits1, vHits1, wHits1, uHits2, vHits2, wHits2;
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const MCParticle *pHitMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            const HitType view{pCaloHit->GetHitType()};
            if (pHitMC == pMC)
            {
                ++primaryHits;
                if (view == HitType::TPC_VIEW_U)
                    uHits1.emplace_back(pCaloHit);
                else if (view == HitType::TPC_VIEW_V)
                    vHits1.emplace_back(pCaloHit);
                else
                    wHits1.emplace_back(pCaloHit);
            }
            else
            {
                ++secondaryHits;
                if (view == HitType::TPC_VIEW_U)
                    uHits2.emplace_back(pCaloHit);
                else if (view == HitType::TPC_VIEW_V)
                    vHits2.emplace_back(pCaloHit);
                else
                    wHits2.emplace_back(pCaloHit);
            }
        }
        std::cout << "  PDG: " << pMC->GetParticleId() << " E: " << pMC->GetEnergy() << " P:" << pMC->GetMomentum().GetMagnitude() <<
            " L: " << (pMC->GetEndpoint() - pMC->GetVertex()).GetMagnitude() << " Primary Hits: " << primaryHits << " Secondary Hits: " <<
            secondaryHits << std::endl;

        std::string suffix{": " + std::to_string(pMC->GetParticleId())};
        if (!uHits1.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits1, "U(1)" + suffix, AUTOITER));
        if (!uHits2.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits2, "U(2)" + suffix, AUTOITER));
        if (!vHits1.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits1, "V(1)" + suffix, AUTOITER));
        if (!vHits2.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits2, "V(2)" + suffix, AUTOITER));
        if (!wHits1.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits1, "W(1)" + suffix, AUTOITER));
        if (!wHits2.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits2, "W(2)" + suffix, AUTOITER));
    }

    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCInfoAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
