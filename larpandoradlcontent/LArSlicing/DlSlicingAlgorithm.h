/**
 *  @file   larpandoracontent/LArMonitoring/DlSlicingAlgorithm.h
 *
 *  @brief  Header file for the pfo validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SLICING_ALGORITHM_H
#define LAR_DL_SLICING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArSliceHelper.h"

namespace lar_content
{

/**
 *  @brief  DlSlicingAlgorithm class
 */
class DlSlicingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlSlicingAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~DlSlicingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Prepare a training sample to train the DL slicing neural network
     *
     *  @return The status code following the attempt to prepare the training sample
     */
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief  Infer using a trained DL slicing neural network
     *
     *  @return The status code following the attempt to infer
     */
    pandora::StatusCode Infer();

    /**
     *  @brief  Populate the root tree with slice information
     *
     *  @param  mcSlices the mc slices
     *  @param  backgroundHits the background hits
     */
    void PopulateRootTree(const LArSliceHelper::SliceHitsMap &mcSlices, const pandora::CaloHitList &backgroundHits) const;

    /**
     *  @brief  Filter the slice hits to only those associated with the cosmic ray MC particle. In the unlikely event that the cosmic ray has no
     *          direct hits, the original slice hits are added to the output list.
     *
     *  @param  sliceHits the slice hits
     *  @param  pCosmicMC the cosmic ray MC particle
     *  @param  cosmicHits the output cosmic hits
     */
    void FilterSliceHitsToCosmic(const pandora::CaloHitList &sliceHits, const pandora::MCParticle *const pCosmicMC, pandora::CaloHitList &cosmicHits) const;

    /**
     *  @brief  Visualize the slices
     *
     *  @param  mcSlices the mc slices
     */
    void VisualizeSlices(const LArSliceHelper::SliceHitsMap &mcSlices) const;

    /**
     *  @brief  Visualize the slices
     *
     *  @param  mcSlices the mc slices
     *  @param  backgroundHits the background hits
     */
    void VisualizeSlices(const LArSliceHelper::SliceHitsMap &mcSlices, const pandora::CaloHitList &backgroundHits) const;

    std::string m_caloHitListName; ///< Name of input calo hit list
    std::string m_rootFileName; ///< Name of output root file
    std::string m_rootTreeName; ///< Name of output root tree
    bool m_trainingMode; ///< Whether to prepare a training sample or infer
    bool m_visualize; ///< Whether to visualize slices
};

} // namespace lar_content

#endif // LAR_DL_SLICING_ALGORITHM_H
