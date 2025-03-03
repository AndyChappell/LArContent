/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_KALMAN_CLUSTER_CREATION_ALGORITHM_H
#define LAR_KALMAN_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArSlicedCaloHitList.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  KalmanClusterCreationAlgorithm class
 */
class KalmanClusterCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    KalmanClusterCreationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~KalmanClusterCreationAlgorithm();

private:
    typedef std::vector<pandora::HitType> ViewVector;
    typedef std::map<pandora::HitType, pandora::OrderedCaloHitList> ViewOrderedHitsMap;
    typedef std::map<const pandora::CaloHit *, pandora::IntVector> HitCandidateCLusterMap;

    /**
     *  @brief  CandidateCluster class
     */
    class CandidateCluster
    {
    public:
        typedef std::tuple<const pandora::CaloHit *, const pandora::CaloHit *, const pandora::CaloHit *> HitTriplet;

        CandidateCluster() = default;

        void AddTriplet(const pandora::CaloHit *const pHit1, const pandora::CaloHit *const pHit2, const pandora::CaloHit *const pHit3);

    private:
        std::vector<HitTriplet> m_hitTriplets;
    };

    /**
     *  @brief  HitAssociation class
     */
    class HitAssociation
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pPrimaryTarget address of the primary target hit
         *  @param  primaryDistanceSquared distance to the primary target hit squared
         */
        HitAssociation(const pandora::CaloHit *const pPrimaryTarget, const float primaryDistanceSquared);

        /**
         *  @brief  Get the primary target
         *
         *  @return the target distance
         */
        const pandora::CaloHit *GetPrimaryTarget() const;

        /**
         *  @brief  Get the primary distance squared
         *
         *  @return the primary distance squared
         */
        float GetPrimaryDistanceSquared() const;

    private:
        const pandora::CaloHit *m_pPrimaryTarget;   ///< the primary target
        float m_primaryDistanceSquared;             ///< the primary distance squared
    };

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Identifies potential clusters. At this stage hits can appear in more than one cluster.
     */
    void IdentifyCandidateClusters(const ViewVector &order);

    /**
     *  @brief  Filter out hits below any MIP threshold and organise into an ordered calo hit list
     *
     *  @param  pCaloHitList input hit list
     *  @param  selectedCaloHitList the output list of selected hits
     */
    pandora::StatusCode FilterCaloHits(const pandora::CaloHitList *const pCaloHitList, pandora::OrderedCaloHitList &selectedCaloHitList) const;

    /**
     *  @brief  Retrieve the event span in x across all views
     *
     *  @param  min the output for the minimum x coordinate
     *  @param  max the output for the maximum x coordinate
     */
    void GetSpanX(float &min, float &max) const;

    /**
     *  @brief  Clean up member variables between events
     */
    pandora::StatusCode Reset() override;

    /**
     *  @brief  Retrieve the slices from a given view/slice map in the vicinty of a given bin. Empty or missing bins will leave the hit lists unfilled
     *
     *  @param  sliceHitMap The slice-to-hit map from which slices should be retrieved
     *  @param  bin The central bin around which slices should be retrieved
     *  @param  caloHits_m1 The output hit list into which hits from the lesser adjacent bin (if populated) will be placed
     *  @param  caloHits_p1 The output hit list into which hits from the greater adjacent bin (if populated) will be placed
     */
    void GetSlices(const LArSlicedCaloHitList::SliceHitMap &sliceHitMap, const size_t bin, pandora::CaloHitVector &caloHits_m1,
        pandora::CaloHitVector &caloHits_p1) const;

    /**
     *  @brief  Retrieve the slices from a given view/slice map in the vicinty of a given bin. Empty or missing bins will leave the hit lists unfilled
     *
     *  @param  sliceHitMap The slice-to-hit map from which slices should be retrieved
     *  @param  bin The central bin around which slices should be retrieved
     *  @param  caloHits_0 The output hit list into which hits from the target bin (if populated) will be placed
     *  @param  caloHits_m1 The output hit list into which hits from the lesser adjacent bin (if populated) will be placed
     *  @param  caloHits_p1 The output hit list into which hits from the greater adjacent bin (if populated) will be placed
     */
    void GetSlices(const LArSlicedCaloHitList::SliceHitMap &sliceHitMap, const size_t bin, pandora::CaloHitVector &caloHits_0,
        pandora::CaloHitVector &caloHits_m1, pandora::CaloHitVector &caloHits_p1) const;

    typedef std::map<pandora::HitType, LArSlicedCaloHitList *> ViewSlicedHitsMap;

    float m_minMipFraction; ///< Minimum fraction of a MIP to consider a hit
    pandora::StringVector m_caloHitListNames; ///< The names of the calo hit lists to cluster
    ViewOrderedHitsMap m_viewHitsMap; ///< Map from the view to the corresponding ordered calo hits
    ViewSlicedHitsMap m_slicedCaloHits; ///< Collection of calo hits in each view organised into slices in the drift coordinate
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline KalmanClusterCreationAlgorithm::HitAssociation::HitAssociation(const pandora::CaloHit *const pPrimaryTarget, const float primaryDistanceSquared) :
    m_pPrimaryTarget(pPrimaryTarget),
    m_primaryDistanceSquared(primaryDistanceSquared)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit *KalmanClusterCreationAlgorithm::HitAssociation::GetPrimaryTarget() const
{
    return m_pPrimaryTarget;
}

inline float KalmanClusterCreationAlgorithm::HitAssociation::GetPrimaryDistanceSquared() const
{
    return m_primaryDistanceSquared;
}

} // namespace lar_content

#endif // #ifndef LAR_KALMAN_CLUSTER_CREATION_ALGORITHM_H
