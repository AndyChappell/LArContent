/**
 *  @file   larpandoracontent/LArWorkshop/MyClusteringAlgorithm.h
 *
 *  @brief  Header file for a custom clustering algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MY_CLUSTERING_ALGORITHM_H
#define LAR_MY_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MyClusteringAlgorithm class
 */
class MyClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MyClusteringAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    size_t m_clusterSize;                              ///< The maximum number of hits to add to a cluster
    pandora::StringVector m_caloHitListNames;       ///< The vector of input calo hit list names
    std::string m_outputListPrefix;                 ///< The prefix to be used for output cluster lists - U/V/W will be added as appropriate
};

} // namespace lar_content

#endif // #ifndef LAR_MY_CLUSTERING_ALGORITHM_H

