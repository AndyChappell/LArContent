/**
 *  @file   larpandoracontent/LArMonitoring/TrackMonitoringAlgorithm.h
 *
 *  @brief  Calculates metrics for protons, charged pions and Michel electrons.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_MONITORING_ALGORITHM_H
#define LAR_TRACK_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  TrackMonitoringAlgorithm class
 */
class TrackMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    TrackMonitoringAlgorithm();

    virtual ~TrackMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Populate the metrics for the MC particles
     *
     *  @param  mcToHitsMap the mapping from MC particles to hits
     *  @param  mcToPfosMap the mapping from MC particles to PFOs
     */
    void PopulateMetrics(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const LArMCParticleHelper::MCToPfosMap &mcToPfosMap) const;

    /**
     *  @brief  Visualize the event, distinguishing hits that are in PFOs and those that are not.
     *
     *  @param  mcToHitsMap the mapping from MC particles to hits
     *  @param  pfoList the list of PFOs in the event
     */
    void Visualize(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const pandora::PfoList &pfoList) const;

    std::string m_filename;         ///< The name of the output ROOT file
    std::string m_treename;         ///< The name of the output ROOT tree
    std::string m_pfoListName;      ///< The name of the PFO list to assess
    std::string m_caloHitListName;  ///< The name of the input calo hit list for all 2D hits
    bool m_visualize;               ///< Whether to visualize events
};

} // namespace lar_content

#endif // LAR_TRACK_MONITORING_ALGORITHM_H
