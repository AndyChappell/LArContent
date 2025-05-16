/**
 *  @file   larpandoracontent/LArMonitoring/TrackOverlapMonitoringAlgorithm.h
 *
 *  @brief  Algorithm for collecting statistics on cases where track overlap near interaction vertices causes reconstruction problems.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_OVERLAP_MONITORING_ALGORITHM_H
#define LAR_TRACK_OVERLAP_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <Eigen/Dense>

namespace lar_content
{

/**
 *  @brief  TrackOverlapMonitoringAlgorithm class
 */
class TrackOverlapMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    TrackOverlapMonitoringAlgorithm();

    virtual ~TrackOverlapMonitoringAlgorithm();

private:
    struct VertexCompare
    {
        bool operator()(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs) const
        {
            const float epsilon{1e-1};

            if (std::fabs(lhs.GetX() - rhs.GetX()) > epsilon)
                return lhs.GetX() < rhs.GetX();
            if (std::fabs(lhs.GetY() - rhs.GetY()) > epsilon)
                return lhs.GetY() < rhs.GetY();
            if (std::fabs(lhs.GetZ() - rhs.GetZ()) > epsilon)
                return lhs.GetZ() < rhs.GetZ();

            return false;
        }
    };

    struct MatrixHit
    {
        bool operator==(const MatrixHit &other) const
        {
            return m_x == other.m_x && m_z == other.m_z;
        }

        float m_x;
        float m_z;
    };

    struct MatrixHitHash
    {
        std::size_t operator()(const MatrixHit &hit) const
        {
            size_t h1 = std::hash<float>()(hit.m_x);
            size_t h2 = std::hash<float>()(hit.m_z);

            return h1 ^ (h2 << 1);
        }
    };

    typedef std::map<const pandora::Pfo *, pandora::MCParticleSet> PfoToMCMap;
    typedef std::map<const pandora::MCParticle *, pandora::PfoSet> MCToPfoMap;
    typedef std::map<const pandora::MCParticle *, pandora::CaloHitSet> MCToHitsMap;
    typedef std::map<const pandora::CartesianVector, pandora::MCParticleSet, VertexCompare> VertexToMCMap;
    typedef std::map<const pandora::MCParticle *, pandora::MCParticleSet> MCToMCMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode CreateMaps();
    pandora::StatusCode FindTrueOverlapCandidates(MCToMCMap &overlapCandidates) const;
    pandora::StatusCode AssessPfos(const MCToMCMap &overlapCandidates, const pandora::PfoList &pfoList) const;
    void CollectHitsByView(const pandora::MCParticle *const pMC, pandora::CaloHitList &uHits, pandora::CaloHitList &vHits, pandora::CaloHitList &wHits) const;
    void VectorizeAndFilterHits(const pandora::CaloHitList &hits, const Eigen::RowVector2f &vertex, const float distance, Eigen::MatrixXf &filteredHits) const;
    void GetDifferenceAndFilterHits(const Eigen::MatrixXf &hits1, const Eigen::MatrixXf &hits2, const float distance, Eigen::MatrixXf &filteredHits1, Eigen::MatrixXf &filteredHits2) const;

    bool m_visualise;   ///< Whether to produce visual monitoring output
    bool m_writeFile;   ///< Whether to produce ROOT output file
    std::string m_rootFileName; ///< The filename of the ROOT output file
    std::string m_rootTreeName; ///< The name of the ROOT tree
    std::string m_pfoListName;  ///< The name of the PFO list to assess
    PfoToMCMap m_pfoToMCMap; ///< PFO to MC map
    MCToPfoMap m_mcToPfoMap; ///< MC to PFO map
    MCToHitsMap m_mcToHitsMap; ///< MC to hits map
    VertexToMCMap m_vertexToMCMap; ///< Vertex to MC map
    float m_vertexRadius; ///< The radius of the vertex to consider
    float m_distance; ///< The distance to consider for overlap
};

} // namespace lar_content

#endif // LAR_TRACK_OVERLAP_MONITORING_ALGORITHM_H
