/**
 *  @file   larpandoracontent/LArTrackShowerId/NnClusterCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the neural network based cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NN_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_NN_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  NnClusterCharacterisationAlgorithm class
 */
class NnClusterCharacterisationAlgorithm : public ClusterCharacterisationBaseAlgorithm
{
public:

private:
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_NN_CLUSTER_CHARACTERISATION_ALGORITHM_H
