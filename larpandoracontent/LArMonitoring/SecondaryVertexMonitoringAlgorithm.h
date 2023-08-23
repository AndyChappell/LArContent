/**
 *  @file   larpandoracontent/LArMonitoring/SecondaryVertexMonitoringAlgorithm.h
 *
 *  @brief  Header file for the secondary vertex monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SECONDARY_VERTEX_MONITORING_ALGORITHM_H
#define LAR_SECONDARY_VERTEX_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SecondaryVertexMonitoringAlgorithm class
 */
class SecondaryVertexMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    SecondaryVertexMonitoringAlgorithm();

    virtual ~SecondaryVertexMonitoringAlgorithm();

private:
    pandora::StatusCode AssessVertices() const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_writeFile;               // Whether to produce ROOT output file
    std::string m_filename;         // The filename of the ROOT output file
    std::string m_treename;         // The name of the ROOT tree
    std::string m_caloHitListName;  // The name of the 2D calo hit list
    std::string m_vertexListName;   // The name of the output secondary vertex list
};

} // namespace lar_content

#endif // LAR_SECONDARY_VERTEX_MONITORING_ALGORITHM_H
