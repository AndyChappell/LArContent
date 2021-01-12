/**
 *  @file   larpandoracontent/LArMonitoring/PfoVisualMonitoringAlgorithm.h
 *
 *  @brief  Header file for the PFO visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_VISUAL_MONITORING_ALGORITHM_H
#define LAR_PFO_VISUAL_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoVisualMonitoringAlgorithm class
 */
class PfoVisualMonitoringAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoVisualMonitoringAlgorithm();

    virtual ~PfoVisualMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_caloHitListName;    ///< Name of input calo hit list
    std::string     m_pfoListName;        ///< Name of input pfo list
    bool            m_foldToPrimaries;    ///< Whether to fold the hierarchy back to primaries
    bool            m_foldToShowers;      ///< Whether to fold the hierarchy back to leading showers
};

} // namespace lar_content

#endif // LAR_PFO_VISUAL_MONITORING_ALGORITHM_H

