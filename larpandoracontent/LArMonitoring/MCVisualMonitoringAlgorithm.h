/**
 *  @file   larpandoracontent/LArMonitoring/MCVisualMonitoringAlgorithm.h
 *
 *  @brief  Header file for the MC visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_VISUAL_MONITORING_ALGORITHM_H
#define LAR_MC_VISUAL_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MCVisualMonitoringAlgorithm class
 */
class MCVisualMonitoringAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCVisualMonitoringAlgorithm();

    virtual ~MCVisualMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_caloHitListName;    ///< Name of input calo hit list
    bool            m_foldToPrimaries;    ///< Whether to fold the hierarchy back to primaries
    bool            m_foldToShowers;      ///< Whether to fold the hierarchy back to leading showers
};

} // namespace lar_content

#endif // LAR_MC_VISUAL_MONITORING_ALGORITHM_H

