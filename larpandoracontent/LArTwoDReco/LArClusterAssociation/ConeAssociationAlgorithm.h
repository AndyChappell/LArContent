/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/ConeAssociationAlgorithm.h
 *
 *  @brief  Header file for the cone association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CONE_ASSOCIATION_ALGORITHM_H
#define LAR_CONE_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ConeAssociationAlgorithm class
 */
class ConeAssociationAlgorithm : public ClusterMergingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ConeAssociationAlgorithm();

    /**
     *  @brief Destructor
     */
    virtual ~ConeAssociationAlgorithm();

private:
    class ViewCluster
    {
    public:
        /**
         *  @brief  Constructs an object that relates clusters from different views that are considered to be candidates for having come
         *          from a common particle.
         *
         *  @param  pClusterU The cluster in the U view
         *  @param  pClusterV The cluster in the V view
         *  @param  pClusterW The cluster in the W view
         *  @param  chi2 The chi-squared value describing the quality of the match between the clusters
         *  @param  overlapStart The low x-cooridnate for the common overlap region
         *  @param  overlapFinish The high x-coordinate for the common overlap region
         */
        ViewCluster(const pandora::Cluster *pClusterU, const pandora::Cluster *pClusterV, const pandora::Cluster *pClusterW,
            const float chi2, const float overlapStart, const float overlapFinish);

        /**
         *  @brief  Retrieves a list of shared views between this ViewCluster and another ViewCluster
         *
         *  @param  other The other ViewCluster to compare against
         *  @param  sharedViews The output list of shared views
         */
        void GetSharedViews(const ViewCluster &other, std::vector<pandora::HitType> &sharedViews) const;

        /**
         *  @brief  Determines if this ViewCluster has at least one shared cluster with another ViewCluster
         *
         *  @param  other The other ViewCluster to compare against
         *  @return true if at least one cluster is common to both ViewClusters, false otherwise
         */
        bool HasSharedCluster(const ViewCluster &other) const;

        /**
         *  @brief  Determines if this ViewCluster shares a cluster in a specified view with another ViewCluster
         *
         *  @param  other The other ViewCluster to compare against
         *  @param  view The view to consider in the comparison
         *  @return true if at least one cluster is common to both ViewClusters, false otherwise
         */
        bool HasSharedCluster(const ViewCluster &other, const pandora::HitType view) const;

        /**
         *  @brief Sort ViewClusters based on size of overlap and chi-squared value.
         *         The primary comparison is between the minimum number of hits from a view in the overlap region, but if these values are
         *         equal, the chi-squared values are compared.
         *
         *  @param  lhs the first object to compare
         *  @param  rhs the second object to compare
         *  @return true if the lhs number of hits is greater than the rhs, false if the rhs is greater, otherwise true if lhs chi-squared
         *          value is lower than rhs, false otherwise
         */
        static bool Sort(const ViewCluster *lhs, const ViewCluster *rhs)
        {
            if (lhs->m_minHitsInOverlapRegion > rhs->m_minHitsInOverlapRegion)
                return true;
            else if (lhs->m_minHitsInOverlapRegion == rhs->m_minHitsInOverlapRegion)
                return lhs->m_chi2 < rhs->m_chi2;
            else
                return false;
        }

        /**
         *  @brief Sort ViewClusters based on size of overlap and chi-squared value.
         *         The primary comparison is between the maximum number of hits from a view in the overlap region, but if these values are
         *         equal, the chi-squared values are compared.
         *
         *  @param  lhs the first object to compare
         *  @param  rhs the second object to compare
         *  @return true if the lhs number of hits is greater than the rhs, false if the rhs is greater, otherwise true if lhs chi-squared
         *          value is lower than rhs, false otherwise
         */
        static bool SortMax(const ViewCluster *lhs, const ViewCluster *rhs)
        {
            if (lhs->m_maxHitsInOverlapRegion > rhs->m_maxHitsInOverlapRegion)
                return true;
            else if (lhs->m_maxHitsInOverlapRegion == rhs->m_maxHitsInOverlapRegion)
                return lhs->m_chi2 < rhs->m_chi2;
            else
                return false;
        }

        const std::string ToString() const
        {
            return "chi2: " + std::to_string(m_chi2) + " hits: " + std::to_string(m_minHitsInOverlapRegion);
        }

        const pandora::ClusterList &GetClusterList() const
        {
            return m_clusters;
        }

        const pandora::Cluster* GetCluster(const pandora::HitType view) const
        {
            for (const pandora::Cluster *pCluster : m_clusters)
                if (LArClusterHelper::GetClusterHitType(pCluster) == view)
                    return pCluster;
            return nullptr;
        }

        float GetChi2() const
        {
            return m_chi2;
        }

    private:
        pandora::ClusterList m_clusters;
        const float m_chi2;
        const float m_overlapStart;
        const float m_overlapFinish;
        const float m_overlapSize;
        int m_minHitsInOverlapRegion;
        int m_maxHitsInOverlapRegion;
    };
    typedef std::map<const pandora::Cluster *, float> ClusterExtremumMap;
    typedef std::map<const pandora::Cluster *, pandora::CaloHitVector> ClusterHitsMap;
    typedef std::vector<const ViewCluster *> ViewClusterVector;
    typedef std::vector<pandora::HitType> HitTypeVector;

    virtual pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Identifies clusters that can be matched across three views
     *
     *  @param  clusterVector the vector of clusters to consider for association
     *  @param  clusterMergeMap the map from seed clusters to candidate clusters for growth
     *  @param  minSeedCaloHits the minimum number of calo hits for a cluster to be considered a seed
     *  @oaram  viewClusterVector The collection of clusters associated across views
     */
    void AssociateClusters(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV,
        const pandora::ClusterList &clusterListW, ViewClusterVector &viewClusterVector) const;

    /**
     *  @brief  Calculate the chi-squared value for the hits in the specified views.
     *
     *  @param  hitsU The hits in the U view
     *  @param  hitsV The hits in the V view
     *  @param  hitsW The hits in the W view
     *  @param  overlapStart The start of the overlap in x
     *  @param  overlapFinish The finish of the overlap in x
     *
     *  @return The chi-squared value for the proposed 3D hits
     */
    float AssessHitOverlap(const pandora::CaloHitVector &hitsU, const pandora::CaloHitVector &hitsV, const pandora::CaloHitVector &hitsW,
        const float overlapStart, const float overlapFinish) const;

    /**
     *  @brief  Adds key cluster details to maps.
     *
     *  @param  pCluster The cluster for which details should be extracted
     *  @param  clusterHitsMap The output map from the cluster to its hits, sorted by x position
     *  @param  clusterXminMap The output map from the cluster to its minimal x coordinate
     *  @param  clusterXmaxMap The output map from the cluster to its maximal x coordinate
     */
    void PopulateClusterMaps(const pandora::Cluster *pCluster, ClusterHitsMap &clusterHitsMap, ClusterExtremumMap &clusterXminMap, 
        ClusterExtremumMap &clusterXmaxMap) const;

    /**
     *  @brief  Takes the set of potential inter-view cluster associations and identifies commonalities between different triplets and
     *          attempts to consolidate them.
     *
     *  @param  viewClusterVector The input vector of ViewClusters
     *  @param  clusterMergeMap The output cluster merge map
     */
    void ConsolidateClusters(ViewClusterVector &viewClusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Takes the set of potential inter-view cluster associations and builds PCA-based cones to determine which other clusters
     *          might be associated with the current triplet.
     *
     *  @param  viewClusterVector The input vector of ViewClusters
     *  @param  clusterMergeMap The output cluster merge map
     */
    void GrowClusters(ViewClusterVector &viewClusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Computes the principal axis of a a cluster.
     *
     *  @param  pCluster The cluster for which the axis is to be computed
     *  @param  p0 The output start coordinate of the axis
     *  @param  p1 The output finish cooridinate of the axis
     *  @param  transverse The transverse axis direction
     *  @param  length The output length of the axis
     *  @param  ratio The output ratio of the 2 principal eigen values to provide a transverse scale relative to the principal axis
     */
    void GetConeParameters(const pandora::Cluster *pCluster, pandora::CartesianVector &start, pandora::CartesianVector &finish,
        pandora::CartesianVector &transverse, float &length, float &ratio) const;

    pandora::StringVector m_inputListNames;     ///< The list of input clusters
    unsigned int m_minClusterLayers;            ///< minimum allowed number of layers for a clean cluster
    mutable int m_runCount;                     ///< the number of times the algorithm has been run
    bool m_visualize;                           ///< whether or not to visualize the algorithm
    std::map<pandora::HitType, std::string> m_viewToListMap;
};

} // namespace lar_content

#endif // #ifndef LAR_CONE_ASSOCIATION_ALGORITHM_H
