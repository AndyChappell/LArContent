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

#include <optional>

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
    typedef std::pair<float, float> ChargeDeposit;
    typedef std::vector<ChargeDeposit> ChargeProfile;
    struct TrajectorySegment
    {
        TrajectorySegment(const pandora::CartesianVector &start, const pandora::CartesianVector &end,
            const ChargeProfile &chargeProfileU, const ChargeProfile &chargeProfileV, const ChargeProfile &chargeProfileW) :
            m_start(start),
            m_end(end),
            m_chargeProfileU(chargeProfileU),
            m_chargeProfileV(chargeProfileV),
            m_chargeProfileW(chargeProfileW)
        {
        }

        pandora::CartesianVector m_start;
        pandora::CartesianVector m_end;
        ChargeProfile m_chargeProfileU;
        ChargeProfile m_chargeProfileV;
        ChargeProfile m_chargeProfileW;
    };
    typedef std::map<int, std::unique_ptr<TrajectorySegment>> Trajectory;

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

    /**
     *  @brief  Estimate the 3D trajectory from the given hits in each view.
     *          Uses the PlaneSolver to identify high-confidence matched 3D points from the 2D hits.
     * 
     *  @param  refHitsU the U view hits
     *  @param  refHitsV the V view hits
     *  @param  refHitsW the W view hits
     *  @param  refPositions3D the resulting ordered collection of 3D points describing the trajectory
     */
    void Estimate3DTrajectory(const pandora::CaloHitList &refHitsU, const pandora::CaloHitList &refHitsV, const pandora::CaloHitList &refHitsW,
        pandora::CartesianPointVector &refPositions3D) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const LArPlaneContextObject *m_planeContextObject;  ///< The PlaneContextObject for the current event
    std::string m_planeSolverContextName;               ///< The name associated with the PlaneContextObject
};

} // namespace lar_content

#endif // #ifndef CALO_MATCH_TOOL_H
