/**
 *  @file   larpandoracontent/LArMonitoring/SliceValidationAlgorithm.h
 *
 *  @brief  Header file for the pfo validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_VALIDATION_ALGORITHM_H
#define LAR_SLICE_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  SliceValidationAlgorithm class
 */
class SliceValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SliceValidationAlgorithm();

private:
    typedef std::map<const pandora::MCParticle *, const pandora::MCParticle *> MCLeadingMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create a map from MCParticles to the calo hits they contribute to
     *
     *  @param  caloHitList the calo hit list
     *  @param  mcToHitsMap the mc to hits map to be filled
     */
    void CreateMCToHitsMap(const pandora::CaloHitList &caloHitList, LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    /**
     *  @brief  Create a map from MCParticles to their leading MCParticle ancestor
     *
     *  @param  mcToHitsMap the mc to hits map
     *  @param  mcToLeadingMap the mc to leading map to be filled
     */
    void CreateMCToLeadingMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCLeadingMap &mcToLeadingMap) const;

    /**
     *  @brief  Create a map from slices to calo hits
     *
     *  @param  mcToHitsMap the mc to hits map
     *  @param  mcToLeadingMap the mc to leading map
     *  @param  sliceToHitsMap the slice to hits map to be filled
     */
    void CreateSliceToHitsMap(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCLeadingMap &mcToLeadingMap,
        LArMCParticleHelper::MCContributionMap &sliceToHitsMap) const;

    std::string m_caloHitListName; ///< Name of input calo hit list
    std::string m_pfoListName; ///< Name of input pfo list
};

} // namespace lar_content

#endif // LAR_SLICE_VALIDATION_ALGORITHM_H
