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

    /**
     *  @brief  Destructor
     */
    ~SliceValidationAlgorithm();

private:
    typedef std::map<const pandora::MCParticle *, const pandora::MCParticle *> MCLeadingMap;

    template <typename Ti, typename Tj>
    using ContingencyTable = std::map<Ti, std::map<Tj, int>>;

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

    /**
     * *  @brief  Validate the slices
     *
     *  @param  mcSlices the mc slices
     *  @param  recoSlices the reconstructed slices
     */
    void ValidateSlices(const LArMCParticleHelper::MCContributionMap &mcSlices, const pandora::PfoList &recoSlices) const;

    typedef std::map<const pandora::MCParticle *, pandora::PfoList> TrueToRecoSliceMap;
    /**
     *  @brief  Match reconstructed slices to true slices
     *
     *  @param  mcSlices the mc slices
     *  @param  recoSlices the reconstructed slices
     *  @param  trueToRecoSliceMap the true to reco slice map to be filled
     */
    void MatchRecoToTrueSlices(const LArMCParticleHelper::MCContributionMap &mcSlices, const pandora::PfoList &recoSlices,
        TrueToRecoSliceMap &trueToRecoSliceMap) const;

    /**
     *  @brief  Build a root tree with various slice reconstruction metrics.
     *          Computes slice efficiency, purity, completeness and ARI, alongside neutrino or cosmic tags and total slice ADC
     *          to permit weighting of metrics by slice size.
     *
     *  @param  trueToRecoSliceMap the true to reco slice map
     *  @param  mcSlices the mc slices
     */
    void PopulateRootTree(const TrueToRecoSliceMap &trueToRecoSliceMap, const LArMCParticleHelper::MCContributionMap &mcSlices) const;

    /**
     *  @brief  Visualize the slices
     *
     *  @param  mcSlices the mc slices
     *  @param  recoSlices the reconstructed slices
     */
    void VisualizeSlices(const LArMCParticleHelper::MCContributionMap &mcSlices, const pandora::PfoList &recoSlices) const;

    /**
     *  @brief  Visualize the slice matches
     *
     *  @param  trueToRecoSliceMap the true to reco slice map
     *  @param  mcSlices the mc slices
     */
    void VisualizeSliceMatches(const TrueToRecoSliceMap &trueToRecoSliceMap, const LArMCParticleHelper::MCContributionMap &mcSlices) const;

    /**
     *  @brief  Find the intersection of two calo hit lists
     *
     *  @param  listA the first calo hit list
     *  @param  listB the second calo hit list
     *  @param  intersection the calo hit list to be filled with the intersection
     */
    void FindSetIntersection(const pandora::CaloHitList &listA, const pandora::CaloHitList &listB,
        pandora::CaloHitList &intersection) const;

    std::string m_caloHitListName; ///< Name of input calo hit list
    std::string m_pfoListName; ///< Name of input pfo list
    std::string m_rootFileName; ///< Name of output root file
    std::string m_rootTreeName; ///< Name of output root tree
};

} // namespace lar_content

#endif // LAR_SLICE_VALIDATION_ALGORITHM_H
