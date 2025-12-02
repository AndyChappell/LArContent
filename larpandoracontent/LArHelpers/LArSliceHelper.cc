/**
 *  @file   larpandoracontent/LArHelpers/LArSliceHelper.cc
 *
 *  @brief  Implementation of the lar monte carlo particle helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArSliceHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

namespace lar_content
{

using namespace pandora;

void LArSliceHelper::GetSliceToHitsMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap,
    const LArMCParticleHelper::MCLeadingMap &mcToLeadingMap, SliceHitsMap &sliceToHitsMap)
{
    for (const auto &[pMC, caloHits] : mcToHitsMap)
    {
        const MCParticle *const pLeadingMC{mcToLeadingMap.at(pMC)};
        // Find the TPCs with hits for this MCParticle
        std::set<unsigned int> tpcIdSet;
        for (const CaloHit *const pCaloHit : caloHits)
        {
            const LArCaloHit *const pLArCaloHit{static_cast<const LArCaloHit *>(pCaloHit)};
            if (pLArCaloHit)
                tpcIdSet.insert(pLArCaloHit->GetLArTPCVolumeId());
        }
        // Initialise the slice entries
        for (const unsigned int tpcId : tpcIdSet)
        {
            const auto key{std::make_pair(tpcId, pLeadingMC)};
            if (sliceToHitsMap.find(key) == sliceToHitsMap.end())
                sliceToHitsMap[key] = CaloHitList();
        }
        // Populate the slice entries
        for (const CaloHit *const pCaloHit : caloHits)
        {
            const LArCaloHit *const pLArCaloHit{static_cast<const LArCaloHit *>(pCaloHit)};
            if (pLArCaloHit)
            {
                const unsigned int tpcId{pLArCaloHit->GetLArTPCVolumeId()};
                sliceToHitsMap[{tpcId, pLeadingMC}].emplace_back(pCaloHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArSliceHelper::FilterSlices(const int pdg, SliceHitsMap &sliceToHitsMap, CaloHitList &backgroundHits)
{
    for (auto it = sliceToHitsMap.begin(); it != sliceToHitsMap.end();)
    {
        const MCParticle *const pMC{it->first.second};
        if (std::abs(pMC->GetParticleId()) == pdg)
        {
            const CaloHitList &hitsToRemove{it->second};
            backgroundHits.insert(backgroundHits.end(), hitsToRemove.begin(), hitsToRemove.end());
            it = sliceToHitsMap.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArSliceHelper::FilterSlices(const IntVector pdgs, SliceHitsMap &sliceToHitsMap, CaloHitList &backgroundHits)
{
    for (const int pdg : pdgs)
        LArSliceHelper::FilterSlices(pdg, sliceToHitsMap, backgroundHits);
}

} // namespace lar_content
