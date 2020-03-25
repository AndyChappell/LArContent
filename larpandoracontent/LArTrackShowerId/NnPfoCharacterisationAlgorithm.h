/**
 *  @file   larpandoracontent/LArTrackShowerId/NnPfoCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the neural network based pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NN_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_NN_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  NnPfoCharacterisationAlgorithm class
 */
class NnPfoCharacterisationAlgorithm : public PfoCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NnPfoCharacterisationAlgorithm();

private:
    bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_NN_PFO_CHARACTERISATION_ALGORITHM_H
