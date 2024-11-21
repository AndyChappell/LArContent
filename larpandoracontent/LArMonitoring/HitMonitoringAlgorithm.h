/**
 *  @file   larpandoracontent/LArMonitoring/HitMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_MONITORING_ALGORITHM_H
#define LAR_HIT_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  HitMonitoringAlgorithm class
 */
class HitMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HitMonitoringAlgorithm();

    virtual ~HitMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the distribution of hit widths in the various views
     *
     *  @param  caloHitList The list of calo hits
     **/
    void GetHitWidthDistribution(const pandora::CaloHitList &caloHitList) const;

    std::string m_filename; ///< The name of the output root file
    std::string m_treename; ///< The name of the output root tree
    std::string m_caloHitListName; ///< The name of the input calo hit list
};

} // namespace lar_content

#endif // LAR_HIT_MONITORING_ALGORITHM_H

