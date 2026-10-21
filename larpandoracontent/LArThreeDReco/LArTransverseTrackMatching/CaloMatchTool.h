/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/CaloMatchTool.h
 *
 *  @brief  Performs matching of clusters between views on the basis of both geometric and calorimetric considerations.
 *          The tool looks for geometric consistency of the hits, but also constructs a normalised profile of charge deposition along the track
 *          in each view, with the profile being the charge deposition along projections of fitted 3D track segments.
 *
 *  $Log: $
 */
#ifndef CALO_MATCH_TOOL_H
#define CALO_MATCH_TOOL_H 1

#include "larpandoracontent/LArObjects/LArPlaneContextObject.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeViewTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CaloMatchTool class
 */
class CaloMatchTool : public TransverseTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CaloMatchTool();

    bool Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    /**
     *  @brief  Construct a 3D trajectory template from the given clusters.
     *          Uses the PlaneSolver to identify high-confidence matched 3D points from the 2D clusters, and build a template trajectory from them.
     * 
     *  @param  pClusterU the U view cluster
     *  @param  pClusterV the V view cluster
     *  @param  pClusterW the W view cluster
     *  @param  trajectory the resulting ordered collection of 3D points describing the trajectory
     *
     *  @return true if the trajectory was successfully constructed, false otherwise
     */
    bool Construct3DTrajectoryTemplates(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
        const pandora::Cluster *const pClusterW, pandora::CartesianPointVector &trajectory) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const LArPlaneContextObject *m_planeContextObject;  ///< The PlaneContextObject for the current event
    std::string m_planeSolverContextName;               ///< The name associated with the PlaneContextObject
};

} // namespace lar_content

#endif // #ifndef CALO_MATCH_TOOL_H
