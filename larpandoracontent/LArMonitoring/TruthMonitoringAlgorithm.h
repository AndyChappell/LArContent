/**
 *  @file   larpandoracontent/LArMonitoring/TruthMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_TRUTH_MONITORING_ALGORITHM_H
#define LAR_TRUTH_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  TruthMonitoringAlgorithm class
 */
class TruthMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    TruthMonitoringAlgorithm();

    virtual ~TruthMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_writeFile;               // Whether to produce ROOT output file
    std::string m_filename;         // The filename of the ROOT output file
    std::string m_treename;         // The name of the ROOT tree
};

} // namespace lar_content

#endif // LAR_TRUTH_MONITORING_ALGORITHM_H
