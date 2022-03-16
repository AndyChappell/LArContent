/**
 *  @file   larpandoracontent/LArUtility/TpcVolumeIteratorAlgorithm.h
 *
 *  @brief  Header file for the TPC volume iterator algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TPC_VOLUME_ITERATOR_ALGORITHM_H
#define LAR_TPC_VOLUME_ITERATOR_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  TpcVolumeIteratorAlgorithm class
 */
class TpcVolumeIteratorAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TpcVolumeIteratorAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_childAlgorithmName;           ///< The name to assign to the child clustering algorithm
    std::string m_caloHitListPrefix;            ///< The prefix for the input CaloHitLists
    std::string m_intermediateClusterListName;  ///< The intermediate name the cild clustering algorithm uses for cluster lists
    std::string m_clusterListPrefix;            ///< The prefix for the output ClusterLists
    pandora::StringVector m_viewVector;         ///< The views to iterate over
    pandora::IntVector m_cryoStatVector;        ///< The cryostats to iterate over
    pandora::IntVector m_tpcVector;             ///< The TPC volumes to iterate over
    pandora::IntVector m_tpcChildVolumeVector;  ///< The TPC child volumes to iterate over
};

} // namespace lar_content

#endif // #ifndef LAR_TPC_VOLUME_ITERATOR_ALGORITHM_H

