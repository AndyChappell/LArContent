/**
 *  @file   larpandoracontent/LArMonitoring/TrackMonitoringAlgorithm.h
 *
 *  @brief  Header file for the track monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_MONITORING_ALGORITHM_H
#define LAR_TRACK_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Match a PFO to MCParticles via their shared hits
     * 
     *  @param  pdg the PDG code of the MCParticle to match to
     *  @param  pfoList the list of PFOs to match
     *  @param  mcToHitsMap map from MCParticles to calo hits
     */
    void MatchPfoToMCParticle(const int pdg, const pandora::PfoList &pfoList, const LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    std::string m_caloHitListName;                       ///< Name of input calo hit list
    std::string m_pfoListName;                           ///< Name of input pfo list
};

} // namespace lar_content

#endif // LAR_TRACK_MONITORING_ALGORITHM_H
