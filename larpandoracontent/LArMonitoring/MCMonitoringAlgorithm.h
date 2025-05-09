/**
 *  @file   larpandoracontent/LArMonitoring/MCMonitoringAlgorithm.h
 *
 *  @brief  Header file for the mc particle monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_MONITORING_ALGORITHM_H
#define LAR_MC_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  MCMonitoringAlgorithm class
 */
class MCMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void PrintHierarchy(const pandora::MCParticle *const pMC, const std::string &tab) const;

    std::string m_caloHitListName;    ///< Name of input calo hit list
    std::string m_mcParticleListName; ///< Name of input MC particle list
};

} // namespace lar_content

#endif // LAR_MC_MONITORING_ALGORITHM_H
