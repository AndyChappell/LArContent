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
     */
    void TagMuonClusters() const;

    /**
     *  @brief  Tag the hits in clusters that aren't muons.
     */
    void TagNonMuonClusters() const;


    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minTrackRatio;                              ///< The minimum ratio of a cluster to the largest in the same MC particle for track tag
    pandora::StringVector m_muonClusterListNames;       ///< The names of the input clusters containing Muons
    pandora::StringVector m_nonMuonClusterListNames;    ///< The names of the input clusters containing particles other than Muons
};

} // namespace lar_content

#endif // #ifndef LAR_HIT_TRUTH_TAGGING_ALGORITHM_H
