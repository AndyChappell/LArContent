/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackFillAlgorithm.h
 *
 *  @brief  Header file for the track recovery algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_FILL_ALGORITHM_H
#define LAR_TRACK_FILL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{

/**
 *  @brief  TrackFillAlgorithm class
 */
class TrackFillAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackFillAlgorithm();

private:
    typedef std::unordered_map<const pandora::CaloHit *, const pandora::Cluster *> HitToClusterToMap;
    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> KDTree;
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> KDNode;
    typedef std::vector<KDNode> NodeList;
    typedef std::unordered_map<const pandora::HitType, KDTree> ViewToKDTreeMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Loops over all clusters to get a map from hits to their existing cluster
     *
     *  @param  clusterList the list of clusters to loop over
     *  @param  hitToClusterMap the map in which to store the mapping of hits to their existing cluster
     */
    void FillHitToClusterMap(const pandora::ClusterList &clusterList, HitToClusterToMap &hitToClusterMap) const;

    /**
     *  @brief  Performs a sliding linear fit to get an ordered set of hits along a cluster trajectory.
     *
     *  @param  pCluster the cluster for which hits should be ordered
     *  @param  orderedHits the list in which to store the ordered hits
     *
     *  @return Whether or not the fit was successful. If false, the ordered hits output is filled with the standard cluster hit order.
     */
    bool FitAndOrderCluster(const pandora::Cluster *const pCluster, pandora::CaloHitList &orderedHits) const;

    /**
     *  @brief  Loops over track-like PFOs and looks for gaps to fill in the tracks. Hits are deemed to fall in the gap if they intercept paths
     *          between hits in the existing clusters. In this case, hits can be moved from their existing PFO to the current PFO.
     *
     *  @param  pfoList the list of track-like PFOs to loop over
     *  @param[in,out]  hitToClusterMap the map of hits to their existing cluster
     *  @param  kdTreeMap the map of KD trees for fast, position-based hit lookup, mapped by view
     */
    void FillTracks(const pandora::PfoList &pfoList, HitToClusterToMap &hitToClusterMap, ViewToKDTreeMap &kdTreeMap) const;

    /**
     *  @brief  Fill KD trees for fast, position-based hit lookup
     *
     *  @param  caloHitList the list of all 2D hits with which to fill the KD trees
     *  @param  kdTreeMap the map in which to store the KD trees, mapped by view
     */
    void FillKDTrees(const pandora::CaloHitList &caloHitList, ViewToKDTreeMap &kdTreeMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitListName; ///< The name of the input list of 2D hits to use for filling the KD trees
    std::string m_inputPfoListName; ///< The name of the input PFO list containing track-like PFOs
    pandora::StringVector m_inputClusterListNames; ///< The names of input cluster lists
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_FILL_ALGORITHM_H
