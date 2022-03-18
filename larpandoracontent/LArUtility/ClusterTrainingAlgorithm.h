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

    /**
     *  @brief  Create child volume hit list
     *
     *  @param  caloHitList  The CaloHitList to persist
     *  @param  listName The name to use for the list
     */
    pandora::StatusCode CreateHitList(const pandora::CaloHitList &caloHitList, std::string listName) const;

    /**
     *  @brief  Write a feature vector for a cluster pairing
     *
     *  @param  view1 The view from which the first set of hits is derived
     *  @param  points1 The set of points from the first cluster
     *  @param  view2 The view from which the second set of hits is derived
     *  @param  points2 The set of points from the second cluster
     *  @param  isSignal Whether or not the feature vector represents signal
     */
    pandora::StatusCode WriteFeatureVector(const pandora::HitType view1, const pandora::CartesianPointVector &points1,
        const pandora::HitType view2, const pandora::CartesianPointVector &points2, const bool isSignal) const;

    bool m_initialisation;          ///< Whether or not this is an initialisation run
    std::string m_trainingFilename; ///< Training file name
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_TRAINING_ALGORITHM_H

