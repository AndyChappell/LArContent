/**
 *  @file   larpandoracontent/LArMonitoring/VertexAssociatedHitMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_ASSOCIATED_HIT_MONITORING_ALGORITHM_H
#define LAR_VERTEX_ASSOCIATED_HIT_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexAssociatedHitMonitoringAlgorithm class
 */
class VertexAssociatedHitMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    VertexAssociatedHitMonitoringAlgorithm();

    virtual ~VertexAssociatedHitMonitoringAlgorithm() = default;

private:
    pandora::StatusCode Run();
    void IdentifyTrackStubs(const pandora::CaloHitList &caloHitList, const pandora::Vertex &vertex) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< The name of the calo hit list
    std::string m_vertexListName;   ///< The name of the vertex list
    float m_transparencyThresholdE; ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float m_energyScaleThresholdE;  ///< Cell energy for which color is at top end of continous color palette
    float m_scalingFactor;          ///< TEve works with [cm], Pandora usually works with [mm] (but LArContent went with cm too)
};

} // namespace lar_content

#endif // LAR_VERTEX_ASSOCIATED_HIT_MONITORING_ALGORITHM_H
