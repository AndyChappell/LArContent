/**
 *  @file   larpandoracontent/LArTrackShowerId/HitTruthTaggingAlgorithm.h
 *
 *  @brief  Header file for the branch growing algorithm base class.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_TRUTH_TAGGING_ALGORITHM_H
#define LAR_HIT_TRUTH_TAGGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  HitTruthTaggingAlgorithm class
 */
class HitTruthTaggingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HitTruthTaggingAlgorithm();

private:
    /**
     *  @brief  Tag the hits in clusters that are muons according to MC truth (in MCC11, delta rays are folded back to the parent muon, we
     *          want to try to separate them).
     *
     *  @param  trackHitList The output list of track-like hits
     *  @param  showerHitList The output list of shower-like hits
     *  @param  diffuseHitList The output list of diffuse-like hits
     */
    void TagMuonClusters(pandora::CaloHitList &trackHitList, pandora::CaloHitList &showerHitList, pandora::CaloHitList &diffuseHitList) const;

    /**
     *  @brief  Tag the hits in clusters that aren't muons.
     *
     *  @param  trackHitList The output list of track-like hits
     *  @param  showerHitList The output list of shower-like hits
     *  @param  michelHitList The output list of Michel-like hits
     *  @param  diffuseHitList The output list of diffuse-like hits
     */
    void TagNonMuonClusters(pandora::CaloHitList &trackHitList, pandora::CaloHitList &showerHitList, pandora::CaloHitList &michelHitList,
        pandora::CaloHitList &diffuseHitList) const;

    /**
     *  @brief  Output the featuer vectors for track/shower ID training.
     *
     *  @param  trackHitList The list of track-like hits
     *  @param  showerHitList The list of shower-like hits
     *  @param  michelHitList The list of Michel-like hits
     *  @param  diffuseHitList The list of diffuse-like hits
     */
    void OutputTruthTagging(const pandora::CaloHitList &trackHitList, const pandora::CaloHitList &showerHitList,
        const pandora::CaloHitList &michelHitList, const pandora::CaloHitList &diffuseHitList) const;

    /**
     *  @brief  Visualise the truth tags for hits.
     *
     *  @param  trackHitList The list of track-like hits
     *  @param  showerHitList The list of shower-like hits
     *  @param  michelHitList The list of Michel-like hits
     *  @param  diffuseHitList The list of diffuse-like hits
     */
    void VisualizeTruthTagging(const pandora::CaloHitList &trackHitList, const pandora::CaloHitList &showerHitList,
        const pandora::CaloHitList &michelHitList, const pandora::CaloHitList &diffuseHitList) const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minTrackRatio;                              ///< The minimum ratio of a cluster to the largest in the same MC particle for track tag
    pandora::StringVector m_muonClusterListNames;       ///< The names of the input clusters containing Muons
    pandora::StringVector m_nonMuonClusterListNames;    ///< The names of the input clusters containing particles other than Muons
    std::string m_outputFileName;                       ///< The base output filename for feature vectors
    bool m_visualize;                                   ///< Whether or not to visualise the truth tagging
    bool m_writeFeatures;                               ///< Whether or not to write the features to training files

    enum Tag : int;
};

} // namespace lar_content

#endif // #ifndef LAR_HIT_TRUTH_TAGGING_ALGORITHM_H
