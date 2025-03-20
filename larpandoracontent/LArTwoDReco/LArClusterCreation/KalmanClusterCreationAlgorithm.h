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
#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include <atomic>
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

    struct KalmanFit
    {
        void InsertHit(const pandora::CaloHit *const pCaloHit);

        const pandora::CaloHit *m_pSeedHit1;
        const pandora::CaloHit *m_pSeedHit2;
        KalmanFilter2D m_kalmanFilter;
        pandora::CaloHitSet m_caloHits;
        const pandora::CaloHit *m_pLastHit;
        int m_id;

        KalmanFit(const pandora::CaloHit *const pSeedHit1, const pandora::CaloHit *const pSeedHit2, const KalmanFilter2D &kalmanFilter, const pandora::CaloHitSet &caloHits,
            const pandora::CaloHit *const pLastHit) :
            m_pSeedHit1(pSeedHit1),
            m_pSeedHit2(pSeedHit2),
            m_kalmanFilter(kalmanFilter),
            m_caloHits(caloHits),
            m_pLastHit(pLastHit),
            m_id(m_counter.fetch_add(1))
        {
        };

    private:
        static std::atomic<int> m_counter;
    };

    typedef std::set<int> KalmanFitIDSet;
    typedef std::vector<KalmanFit> KalmanFitVector;
    typedef std::map<const pandora::CaloHit *, KalmanFitIDSet> HitKalmanFitMap;

    /**
     *  @brief  CandidateCluster class
     */
    class CandidateCluster
    {
    public:
        typedef std::tuple<const pandora::CaloHit *, const pandora::CaloHit *, const pandora::CaloHit *> HitTriplet;
        typedef std::vector<HitTriplet> HitTripletVector;

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
     *
     *  @param  order the orderied vector of views
     */
    void IdentifyCandidateClusters(const ViewVector &order);

    /**
     *  @brief  Make cluster seeds from hits in a single view
     *
     *  @param[in]  sliceCaloHits the input slice calo hits
     *  @param[out]  kalmanFits the vector of kalman fits to update
     *  @param[out]  hitKalmanFitMap the map from hits to kalman fits to update
     */
    void MakeClusterSeeds(const pandora::CaloHitVector &sliceCaloHits, KalmanFitVector &kalmanFits, HitKalmanFitMap &hitKalmanFitMap);

    /**
     *  @brief  Build clusters from hits in a single view
     *
     *  @param[in]  sliceCaloHits the input slice calo hits
     *  @param[in,out]  kalmanFits the vector of kalman fits to update
     *  @param[in,out]  hitKalmanFitMap the map from hits to kalman fits to update
     */
    void BuildClusters(const pandora::CaloHitVector &sliceCaloHits, KalmanFitVector &kalmanFits, HitKalmanFitMap &hitKalmanFitMap);

    /**
     *  @brief  Remove duplicate kalman fits
     *
     *  @param[in,out]  kalmanFits the vector of kalman fits to update
     *  @param[in,out]  hitKalmanFitMap the map from hits to kalman fits to update
     */
    void RemoveDuplicateKalmanFits(KalmanFitVector &kalmanFits, HitKalmanFitMap &hitKalmanFitMap);

    /**
     *  @brief  Make candidate 3D hits from hits in two or more views
     *
     *  @param  caloHits0 the first vector of input 2D calo hits
     *  @param  caloHits1 the second vector of input 2D calo hits
     *  @param  caloHits2 the third vector of input 2D calo hits
     *  @param  triplets the output triplets (or doublets with element 3 null)
     */
    void Make3DHitPermutations(const pandora::CaloHitVector &caloHits0, const pandora::CaloHitVector &caloHits1, const pandora::CaloHitVector &caloHits2, CandidateCluster::HitTripletVector &triplets) const;

    /**
     *  @brief  Filters a collection of 3D hit triplets (or doublets) based on the quality of the match, as determined by the chi-squared value
     *
     *  @param  triplets the input/output triplets (or doublets with element 3 null)
     *  @param  hits3D the output 3D hit positions associated with the 2D hits
     *  @param  chi2s the output chi-squared values for the 3D hita
     */
    void Filter3DHitPermutations(CandidateCluster::HitTripletVector &triplets, pandora::CartesianPointVector &hits3D, pandora::FloatVector &chi2s) const;

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
     *  @param  caloHits The output hit list into which hits from the target bin and adjacent bins (if populated) will be placed
     */
    void GetSlices(const LArSlicedCaloHitList::SliceHitMap &sliceHitMap, const size_t bin, pandora::CaloHitVector &caloHits) const;

    /**
     *  @brief  Determines if two hits are within a given proximity of each other
     *
     *  @param  pCaloHit1 the first calo hit
     *  @param  pCaloHit2 the second calo hit
     *  @param  proximity the radius of the neighbourhood considered proximate
     *
     *  @return true if the hits are within a given proximity of each other, false otherwise
     */
    bool Proximate(const pandora::CaloHit *const pCaloHit1, const pandora::CaloHit *const pCaloHit2, const float proximity = 1.f) const;

    /**
     *  @brief  Determines if a position is contained within a hit
     *
     *  @param  pCaloHit the calo hit
     *  @param  position the position to check
     *  @param  xTol the containment tolerance in the x direction
     *  @param  zTol the containment tolerance in the z direction
     *
     *  @return true if the position is contained within the hit, false otherwise
     */
    bool Contains(const pandora::CaloHit *const pCaloHit, const Eigen::VectorXd &position, const float xTol = 0.f, const float zTol = 0.f) const;

    /**
     *  @brief  Determines if there is a hit between two hits. If the intervening hit significantly overlaps with either of the extremal hits
     *          then the two hits are not considered to skip over the intervening hit.
     *
     *  @param  caloHits the collection of hits within the slice
     *  @param  pCaloHit1 the first calo hit
     *  @param  pCaloHit2 the second calo hit
     *
     *  @return true if there is a hit between the two hits, false otherwise
     */
    bool SkipsOverHit(const pandora::CaloHitVector &caloHits, const pandora::CaloHit *const pCaloHit1, const pandora::CaloHit *const pCaloHit2) const;

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
