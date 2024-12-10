/**
 *  @file   larpandoracontent/LArMonitoring/CheatingClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_CHEATING_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingClusterSplittingAlgorithm class
 */
class CheatingClusterSplittingAlgorithm : public pandora::Algorithm
{
public:
    /**
    *  @brief  Default constructor
    */
    CheatingClusterSplittingAlgorithm();

    virtual ~CheatingClusterSplittingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_visualize;                ///< Whether to produce visual monitoring output
    std::string m_clusterListName;   ///< The name of the input cluster list
    std::string m_caloHitListName;   ///< The name of the input 2D calo hit list
};

} // namespace lar_content

#endif // LAR_CHEATING_CLUSTER_SPLITTING_ALGORITHM_H

