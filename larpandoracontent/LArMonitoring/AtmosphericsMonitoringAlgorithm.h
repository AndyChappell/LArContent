/**
 *  @file   larpandoracontent/LArMonitoring/AtmosphericsMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_ATMOSPHERICS_MONITORING_ALGORITHM_H
#define LAR_ATMOSPHERICS_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  AtmosphericsMonitoringAlgorithm class
 */
class AtmosphericsMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    AtmosphericsMonitoringAlgorithm();

    virtual ~AtmosphericsMonitoringAlgorithm();

private:
    pandora::StatusCode ConstructRootTree() const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_outputFileName; ///< Output file name for training examples
    std::string m_outputTreeName; ///< Output vertex list name
};

} // namespace lar_content

#endif // LAR_ATMOSPHERICS_MONITORING_ALGORITHM_H
