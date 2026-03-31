/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.h
 *
 *  @brief  Header file for the track recovery algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
#define LAR_TRACK_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <optional>

namespace lar_content
{

/**
 *  @brief  TrackRecoveryAlgorithm class
 */
class TrackRecoveryAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackRecoveryAlgorithm();

private:
    typedef std::unordered_map<pandora::HitType, const pandora::Cluster *> ViewToClusterMap;
    typedef std::unordered_map<const pandora::Pfo *, ViewToClusterMap> PfoToViewClusterMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::CaloHitList> ClusterToHitMap;
    typedef std::unordered_map<const pandora::Cluster *, const TwoDSlidingFitResult> ClusterToFitMap;
    typedef std::unordered_map<const pandora::CaloHit *, const pandora::Cluster *> HitToClusterToMap;
    typedef std::unordered_map<const pandora::Cluster *, const pandora::Pfo *> ClusterToPfoMap;
    typedef std::unordered_map<pandora::HitType, pandora::CaloHitList> ViewToHitsMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Performs a sliding linear fit to get an ordered set of hits along a cluster trajectory.
     *
     *  @param  pCluster the cluster for which hits should be ordered
     *  @param  orderedHits the list in which to store the ordered hits
     */
    std::optional<TwoDSlidingFitResult> FitAndOrderCluster(const pandora::Cluster *const pCluster, pandora::CaloHitList &orderedHits) const;

    /**
     *  @brief  Identifies hits that are missing from the PFO based on independent triplet matching. hitsA acts as the reference list to look
     *          up the hit triplets. Hits are deemed missing if two of the three hits in a triplet are present in the PFO.
     *
     *  @param  hitsA the list of hits in the reference view
     *  @param  hitsB the list of hits in the second view
     *  @param  hitsC the list of hits in the third view
     *  @param  unmatchedHitsB the set in which to store hits missing in only the second view
     *  @param  unmatchedHitsC the set in which to store hits missing in only the third view
     */
    void FindUnmatchedHits(const pandora::CaloHitList &hitsA, const pandora::CaloHitList &hitsB, const pandora::CaloHitList &hitsC,
        pandora::CaloHitSet &unmatchedHitsB, pandora::CaloHitSet &unmatchedHitsC) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputPfoListName; ///< The name of the input PFO list containing track-like PFOs
    pandora::StringVector m_inputClusterListNames; ///< The list of cluster list names to consider when looking for clusters to recover
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_RECOVERY_ALGORITHM_H
