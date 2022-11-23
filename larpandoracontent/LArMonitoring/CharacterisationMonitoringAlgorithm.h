/**
 *  @file   larpandoracontent/LArMonitoring/CharacterisationMonitoringAlgorithm.h
 *
 *  @brief  Header file for the characterisation monitoring algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHARACTERISATION_MONITORING_ALGORITHM_H
#define LAR_CHARACTERISATION_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_content
{

/**
 *  @brief  CharacterisationMonitoringAlgorithm class
 */
class CharacterisationMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CharacterisationMonitoringAlgorithm();

    /**
     *  @brief Destructor
     */
    virtual ~CharacterisationMonitoringAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_minHitsForGoodView; ///< The minimum number of hits for a good view
    std::string m_trackPfoListName; ///< Name of input track PFO list
    std::string m_showerPfoListName; ///< Name of input shower PFO list
    std::string m_rootTreeName;               ///< The ROOT tree name
    std::string m_rootFileName;               ///< The ROOT file name
};

} // namespace lar_content

#endif // #ifndef LAR_CHARACTERISATION_MONITORING_ALGORITHM_H

