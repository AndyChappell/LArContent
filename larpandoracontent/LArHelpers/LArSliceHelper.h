/**
 *  @file   larpandoracontent/LArHelpers/LArSliceHelper.h
 *
 *  @brief  Header file for the lar monte carlo particle helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_HELPER_H
#define LAR_SLICE_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  LArSliceHelper class
 */
class LArSliceHelper
{
public:
    typedef std::pair<unsigned int, const pandora::MCParticle *> SliceIdentifier;
    typedef std::map<SliceIdentifier, pandora::CaloHitList> SliceHitsMap;

    /**
     *  @brief  Create a map from slices to calo hits
     *
     *  @param  mcToHitsMap the mc to hits map
     *  @param  mcToLeadingMap the mc to leading map
     *  @param  sliceToHitsMap the slice to hits map to be filled
     */
    static void GetSliceToHitsMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const LArMCParticleHelper::MCLeadingMap &mcToLeadingMap,
        SliceHitsMap &sliceToHitsMap);
};

} // namespace lar_content

#endif // #ifndef LAR_SLICE_HELPER_H
