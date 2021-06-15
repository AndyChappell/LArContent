/**
 *  @file   larpandoracontent/LArMonitoring/ShowerCaloMonitoringAlgorithm.h
 *
 *  @brief  Header file for the shower calo algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_CALO_MONITORING_ALGORITHM_H
#define LAR_SHOWER_CALO_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  ShowerCaloMonitoringAlgorithm class
 */
class ShowerCaloMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerCaloMonitoringAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~ShowerCaloMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MakeSelection(const pandora::MCParticleList *pMCList, const pandora::CaloHitList *pCaloHitList,
        LArMCParticleHelper::MCContributionMap &mcMap) const;

    std::string m_caloHitListName; ///< Name of input calo hit list
    bool m_visualize;              ///< Whether or not to visualize MC particles
};

} // namespace lar_content

#endif // LAR_SHOWER_CALO_MONITORING_ALGORITHM_H
