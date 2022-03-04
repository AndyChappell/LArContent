/**
 *  @file   larpandoracontent/LArUtility/ClusterTrainingAlgorithm.h
 *
 *  @brief  Header file for the cluster training algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_TRAINING_ALGORITHM_H
#define LAR_CLUSTER_TRAINING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterTrainingAlgorithm class
 */
class ClusterTrainingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterTrainingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_TRAINING_ALGORITHM_H

