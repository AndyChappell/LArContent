/**
 *  @file   larpandoracontent/LArMonitoring/NeutrinoMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_MONITORING_ALGORITHM_H
#define LAR_NEUTRINO_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoMonitoringAlgorithm class
 */
class NeutrinoMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoMonitoringAlgorithm();

    virtual ~NeutrinoMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Selects the MC particles to consider based on reconstructability criteria
     *
     *  @param  calotHitList The input calo hit list
     *  @param  mcMap The output map from MC particles to calo hits
     **/
    void MakeSelection(const pandora::CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcMap) const;

    std::string m_caloHitListName;  ///< Name of input calo hit list
};

} // namespace lar_content

#endif // LAR_NEUTRINO_MONITORING_ALGORITHM_H
