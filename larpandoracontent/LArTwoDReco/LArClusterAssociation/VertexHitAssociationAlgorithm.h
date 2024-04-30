/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/VertexHitAssociationAlgorithm.h
 *
 *  @brief  Header file for the longitudinal association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_HIT_ASSOCIATION_ALGORITHM_H
#define LAR_VERTEX_HIT_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Helpers/ClusterFitHelper.h"

namespace lar_content
{

/**
 *  @brief  VertexHitAssociationAlgorithm class
 */
class VertexHitAssociationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexHitAssociationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;       ///< The name of the list containing unassociated hits
    unsigned int m_minClusterLayers;     ///< minimum allowed number of layers for a clean cluster
    unsigned int m_maxGapLayers;         ///< maximum allowed number of layers between associated clusters
    unsigned int m_fitLayers;            ///< number of layers to fit at start and end of cluster
    float m_maxGapDistanceSquared;       ///< maximum allowed distance (squared) between associated clusters
    float m_minCosRelativeAngle;         ///< maximum allowed relative angle between associated clusters
    float m_maxTransverseDisplacement;   ///< maximum allowed transverse displacement after extrapolation (normalised to cell size)
    float m_maxLongitudinalDisplacement; ///< maximum allowed longitudinal displacement after extrapolation (normalised to cell size)
    float m_hitSizeZ;                    ///< estimated hit size in z (wire number) dimension, units cm
    float m_hitSizeX;                    ///< estimated hit size in x (drift time) dimension, units cm
    mutable pandora::HitType m_view;     ///< The view to which the hits under consideration belong
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_HIT_ASSOCIATION_ALGORITHM_H
