/**
 *  @file   larpandoracontent/LArMonitoring/PiZeroMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_PIZERO_MONITORING_ALGORITHM_H
#define LAR_PIZERO_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraEnumeratedTypes.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

namespace lar_content
{

/**
 *  @brief  PiZeroMonitoringAlgorithm class
 */
class PiZeroMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PiZeroMonitoringAlgorithm();

    virtual ~PiZeroMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< Name of input calo hit list
    std::string m_pfoListName;      ///< Name of input PFO list
};

} // namespace lar_content

#endif // LAR_PIZERO_MONITORING_ALGORITHM_H
