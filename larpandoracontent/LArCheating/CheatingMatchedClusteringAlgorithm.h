/**
 *  @file   larpandoracontent/LArMonitoring/CheatingMatchedClusteringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_MATCHED_CLUSTERING_ALGORITHM_H
#define LAR_CHEATING_MATCHED_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingMatchedClusteringAlgorithm class
 */
class CheatingMatchedClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
    *  @brief  Default constructor
    */
    CheatingMatchedClusteringAlgorithm();

    virtual ~CheatingMatchedClusteringAlgorithm();

private:
    typedef std::vector<std::string> StringVector;
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_visualize; ///< Whether to produce visual monitoring output
    std::string m_seedClusterListName; ///< The name of the list containing the seed clusters
    StringVector m_caloHitListNames; ///< The names of the input calo hit lists that need to be matched
};

} // namespace lar_content

#endif // LAR_CHEATING_MATCHED_CLUSTERING_ALGORITHM_H

