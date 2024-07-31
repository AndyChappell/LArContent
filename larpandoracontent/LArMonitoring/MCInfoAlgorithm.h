/**
 *  @file   larpandoracontent/LArMonitoring/MCInfoAlgorithm.h
 *
 *  @brief  Header file for the mc particle monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_INFO_ALGORITHM_H
#define LAR_MC_INFO_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  MCInfoAlgorithm class
 */
class MCInfoAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCInfoAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;    ///< Name of input calo hit list
};

} // namespace lar_content

#endif // LAR_MC_INFO_ALGORITHM_H
