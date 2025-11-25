/**
 *  @file   larpandoracontent/LArMonitoring/SliceValidationAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/SliceValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SliceValidationAlgorithm::SliceValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceValidationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    this->CreateMCToHitsMap(*pCaloHitList, mcToHitsMap);

    MCLeadingMap mcToLeadingMap;
    this->CreateMCToLeadingMap(mcToHitsMap, mcToLeadingMap);

    LArMCParticleHelper::MCContributionMap sliceToHitsMap;
    this->CreateSliceToHitsMap(mcToHitsMap, mcToLeadingMap, sliceToHitsMap);
    
    for (const auto &[pSliceMC, caloHits] : sliceToHitsMap)
    {
        std::cout << "Slice MC Particle: PDG=" << std::abs(pSliceMC->GetParticleId()) << " E=" << pSliceMC->GetEnergy()
                  << " NHits=" << caloHits.size() << std::endl;
    }

    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::CreateMCToHitsMap(const CaloHitList &caloHitList, LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            if (mcToHitsMap.find(pMC) == mcToHitsMap.end())
                mcToHitsMap[pMC] = CaloHitList();
            mcToHitsMap[pMC].emplace_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::CreateMCToLeadingMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCLeadingMap &mcToLeadingMap) const
{
    for (const auto &[pMC, _] : mcToHitsMap)
    {
        const MCParticle *pParent{pMC};
        while (!pParent->GetParentList().empty())
            pParent = pParent->GetParentList().front();
        mcToLeadingMap[pMC] = pParent;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceValidationAlgorithm::CreateSliceToHitsMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCLeadingMap &mcToLeadingMap,
    LArMCParticleHelper::MCContributionMap &sliceToHitsMap) const
{
    for (const auto &[pMC, caloHits] : mcToHitsMap)
    {
        const MCParticle *const pLeadingMC{mcToLeadingMap.at(pMC)};
        if (sliceToHitsMap.find(pLeadingMC) == sliceToHitsMap.end())
            sliceToHitsMap[pLeadingMC] = CaloHitList();
        sliceToHitsMap[pLeadingMC].insert(sliceToHitsMap.at(pLeadingMC).end(), caloHits.begin(), caloHits.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
