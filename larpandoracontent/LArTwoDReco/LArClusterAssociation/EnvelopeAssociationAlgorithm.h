/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/EnvelopeAssociationAlgorithm.h
 *
 *  @brief  Header file for the proximity association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H
#define LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EnvelopeAssociationAlgorithm class
 */
class EnvelopeAssociationAlgorithm : public ClusterMergingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EnvelopeAssociationAlgorithm();

private:
    struct Cone
    {
        /**
         *  @brief Constructor
         *
         *  @param  a a vertex describing the cone
         *  @param  b a vertex describing the cone
         *  @param  c a vertex describing the cone
         */
        Cone(const pandora::CartesianVector &a, const pandora::CartesianVector &b, const pandora::CartesianVector &c) :
            m_a(a),
            m_b(b),
            m_c(c)
        {
        }

        const pandora::CartesianVector m_a; ///< A vertex of the bounding envelope
        const pandora::CartesianVector m_b; ///< A vertex of the bounding envelope
        const pandora::CartesianVector m_c; ///< A vertex of the bounding envelope
    };

    class AssociationCandidate
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  cone The bounding region for hits
         *  @param  clusters the clusters contained within the bounding region
         */
        AssociationCandidate(const Cone &cone, const pandora::ClusterList &clusters);

        /**
         *  @brief  Retrieve the list of associated clusters
         *
         *  @return The list of assoiciated clusters
         */
        const pandora::ClusterList &GetClusters() const
        {
            return m_clusters;
        }

        /**
         *  @brief  Retrieve the hits for this association candidate
         *
         *  @return The hits for the association candidate
         */
        const pandora::CaloHitList &GetCaloHits() const
        {
            return m_hits;
        }

        /**
         *  @brief  Retrieve the sorted hits for this association candidate
         *
         *  @return The sorted hits for the association candidate
         */
        const pandora::CaloHitVector &GetSortedCaloHits() const
        {
            return m_sortedHits;
        }

        /**
         *  @brief  Retrieve the minimum x coordinate
         *
         *  @return The minimum x coordinate for the candidate
         */
        float GetMinX() const
        {
            return m_xMin;
        }

        /**
         *  @brief  Retrieve the maximum x coordinate
         *
         *  @return The maximum x coordinate for the candidate
         */
        float GetMaxX() const
        {
            return m_xMax;
        }

        /**
         *  @brief  Retrieve the minimum and maximum overlap fractions between this association candidate and another.
         *          A candidate whose x-span is completely contained within the middle one-third of the x-span of another candidate would
         *          have an overlap fraction of 1 (i.e. it is completely contained), but the other candidate would see it shares only one
         *          third of its span with the other candidate and thus have an overlap fraction of 0.33.
         *
         *  @param  other the candidate against which overlap should be calculated
         *  @param  minOverlapFraction the output minimum overlap fraction
         *  @param  maxOverlapFraction the output maximum overlap fraction
         */
        void GetOverlapFraction(const AssociationCandidate &other, float &minOverlapFraction, float &maxOverlapFraction) const;

    private:
        Cone m_cone;                            ///< The bounding envelope
        pandora::ClusterList m_clusters;        ///< The clusters contained within the bounding envelope
        pandora::CaloHitList m_hits;            ///< The hits belonging to the contained clusters
        pandora::CaloHitVector m_sortedHits;    ///< The sorted hits belonging to the contained clusters
        float m_xMin;                           ///< The minimum x coordinate of the contained hits
        float m_xMax;                           ///< The maximum x coordinate of the contained hits
    };
    typedef std::list<AssociationCandidate> AssociationCandidateList;
    typedef std::map<pandora::HitType, AssociationCandidateList> ViewToAssociationCandidatesMap;
    typedef std::map<pandora::HitType, AssociationCandidate> ViewToAssociationCandiateMap;

    class OverlapCandidates
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pandora the Pandora instance
         *  @param  candidate1 a candidate set of clusters for merging in one view
         *  @param  candidate2 a candidate set of clusters for merging in a second view
         */
        OverlapCandidates(const pandora::Pandora &pandora, const AssociationCandidate &candidate1, const AssociationCandidate &candidate2);

        /**
         *  @brief  Constructor
         *
         *  @param  pandora the Pandora instance
         *  @param  candidate1 a candidate set of clusters for merging in one view
         *  @param  candidate2 a candidate set of clusters for merging in a second view
         *  @param  candidate3 a candidate set of clusters for merging in a third view
         */
        OverlapCandidates(const pandora::Pandora &pandora, const AssociationCandidate &candidate1, const AssociationCandidate &candidate2,
            const AssociationCandidate &candidate3);

        /**
         *  @brief Adds these overlap candidates to the merge map
         *
         *  @param  clusterMergeMap the cluster merge map to populate
         */
        void AddToMergeMap(ClusterMergeMap &clusterMergeMap) const;

        /**
         *  @brief Check if this object shares a seed cluster with another OverlapCandidates object
         *
         *  @param other the other OverlapCandidates object
         */
        bool SharesSeed(const OverlapCandidates &other) const;

        /**
         *  @brief Retrieve the chi-squared value describing the quality of the match between the two candidates.
         *
         *  @return the chi-squared value
         */
        float GetChiSquared() const
        {
            return m_chi2;
        }

        // ATTN: Current 3 v 2 view handling is too favourable to the three view case, need to understand the 3 v 2 view implications and apply
        // some kind of weighting that favours the 3 view case when things are close, but not if the 3 view match is junk
        /**
         *  @brief Check if this object has a lower chi-squared value than the other object. If one candidate has more views than the other,
         *         favour that candidate.
         *
         *  @param  other the other object against which to compare
         *  @return true if the chi-squared value is lower, false otherwise
         */
        bool operator <(const OverlapCandidates &other) const
        {
            if (m_nViews == other.m_nViews)
                return m_chi2 > other.m_chi2;
            else
                return m_nViews > other.m_nViews;
        }

    private:
        int m_nViews;                               ///< The number of views considered
        pandora::ClusterList m_candidateClusters1;  ///< An association candidate in one view
        pandora::ClusterList m_candidateClusters2;  ///< An association candidate in a second view
        pandora::ClusterList m_candidateClusters3;  ///< An association candidate in a third view
        float m_chi2;                               ///< The chi-squared value for the match between views
    };

    virtual pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Identifies good seed clusters and then attempts to grow them based on PCA-derived bounding regions
     *
     *  @param  clusterVector the vector of clusters to consider for association
     *  @param  clusterMergeMap the map from seed clusters to candidate clusters for growth
     *  @param  minSeedCaloHits the minimum number of calo hits for a cluster to be considered a seed
     */
    void AssociateClusters(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap, const unsigned int minSeedCaloHits) const;

    /**
     *  @brief  Determine whether a target cluster is contained by a bounding region
     *
     *  @param  cone the bounding cone
     *  @param  pCluster the target cluster
     *
     *  @return whether the cluster is contained within the bounding region
     */
    bool IsClusterContained(const Cone &cone, const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Determines if a point is contained within a triangle using the Barycentric technique, as described in
     *          Real-Time Collision Detection, C. Ericson (Morgan Kaufmann, 2005).
     *
     *  @param  ab the vector describing the triangle edge ab
     *  @param  ac the vector describing the traingle edge ac
     *  @param  ap the vector describing the displacement of the point p from the point a
     *
     *  @return true if the point is contained within the triangle, false otherwise
     */
    bool IsPointContained(const pandora::CartesianVector &ab, const pandora::CartesianVector &ac, const pandora::CartesianVector &ap) const;

    /**
     *  @brief  Retrieve the cone the bounding envelope to be used during cluster association
     *
     *  @param  pCluster the cluster for which the bounding envelope should be constructed
     *  @return the cone describing the bounding region
     */
    Cone GetBoundingCone(const pandora::Cluster *const pCluster) const;

    pandora::StringVector m_inputListNames;     ///< The list of input clusters
    unsigned int m_minClusterLayers;            ///< minimum allowed number of layers for a clean cluster
    unsigned int m_minSeedCaloHits;             ///< minimum number of calo hits to form a seed for PCA
    mutable int m_runCount;                     ///< the number of times the algorithm has been run
    bool m_visualize;                           ///< whether or not to visualize the algorithm
    std::map<pandora::HitType, std::string> m_viewToListMap;
};

} // namespace lar_content

#endif // #ifndef LAR_ENVELOPE_ASSOCIATION_ALGORITHM_H
