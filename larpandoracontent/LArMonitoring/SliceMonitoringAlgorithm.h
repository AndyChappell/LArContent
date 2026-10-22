/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_MONITORING_ALGORITHM_H
#define LAR_SLICE_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SliceMonitoringAlgorithm class
 */
class SliceMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    SliceMonitoringAlgorithm();

    virtual ~SliceMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;   ///< The name of the input calo hit list
    std::string m_vertexListName;    ///< The name of the reconstructed neutrino vertex list
};

} // namespace lar_content

#endif // LAR_SLICE_MONITORING_ALGORITHM_H
