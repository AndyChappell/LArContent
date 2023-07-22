/**
 *  @file   larpandoracontent/LArWorkshop/MyParticleCreationAlgorithm.h
 *
 *  @brief  Header file for a custom particle creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MY_PARTICLE_CREATION_ALGORITHM_H
#define LAR_MY_PARTICLE_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  MyParticleCreationAlgorithm class
 */
class MyParticleCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MyParticleCreationAlgorithm();

private:
    typedef std::map<pandora::HitType, pandora::ClusterVector> ClusterMap;

    pandora::StatusCode Run();

    float GetFigureOfMerit(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV,
        const TwoDSlidingFitResult &fitResultW) const;

    pandora::StatusCode GetLongClusters(const std::string &clusterListName, pandora::ClusterVector &sortedClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_fitWindow;      ///< The sliding fit window size
    int m_nSamplingPoints;  ///< The number of sampling points to use in fit comparisons
    pandora::StringVector m_clusterListNames;    ///< The names of the cluster lists to process
    std::string m_pfoListName; ///< The output PFO list name
};

} // namespace lar_content

#endif // #ifndef LAR_MY_PARTICLE_CREATION_ALGORITHM_H

