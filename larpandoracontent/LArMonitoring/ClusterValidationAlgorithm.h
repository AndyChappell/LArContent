/**
 *  @file   larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_VALIDATION_ALGORITHM_H
#define LAR_CLUSTER_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterValidationAlgorithm class
 */
class ClusterValidationAlgorithm : public pandora::Algorithm
{
private:
    struct ClusterMetrics
    {
        ClusterMetrics();

        const pandora::MCParticle *m_pMC;
        int m_nContributions;
        float m_wRecoHits;
        float m_purity;
        float m_completeness;
        float m_fragmentationFraction;
    };

public:
    /**
    *  @brief  Default constructor
    */
    ClusterValidationAlgorithm();

    virtual ~ClusterValidationAlgorithm();

private:
    typedef std::vector<std::string> StringVector;
    typedef std::map<int, pandora::ClusterList> ViewClustersMap;
    typedef std::map<const pandora::Cluster *, ClusterMetrics> ClusterMetricsMap;
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Retrieve the metrics of every cluster
     *
     *  @param viewClusterMap The clusters (mapped by view) to compute the purity for
     *  @param metricsMap The output map from clusters to their associated metrics
     */
    void GetMetrics(const ViewClustersMap &viewClusterMap, ClusterMetricsMap &metricsMap) const;

    bool m_visualize;                ///< Whether to produce visual monitoring output
    bool m_writeFile;                ///< Whether to produce ROOT output file
    std::string m_fileName;          ///< The filename of the ROOT output file
    std::string m_treeName;          ///< The name of the ROOT tree
    std::string m_caloHitListName;   ///< The name of the hit list containing all hits in all views
    StringVector m_clusterListNames; ///< The names of the lists of clusters to process
};

} // namespace lar_content

#endif // LAR_CLUSTER_VALIDATION_ALGORITHM_H

