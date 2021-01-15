/**
 *  @file   larpandoracontent/LArMonitoring/MCVisualMonitoringAlgorithm.h
 *
 *  @brief  Header file for the MC visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_VISUAL_MONITORING_ALGORITHM_H
#define LAR_MC_VISUAL_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  MCVisualMonitoringAlgorithm class
 */
class MCVisualMonitoringAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCVisualMonitoringAlgorithm();

    virtual ~MCVisualMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Visualise the MC particles according to their PDG codes
     *
     *  @param  mcMap The map from MC particles to calo hits
     **/
    void VisualiseByPdgCode(const LArMCParticleHelper::MCContributionMap &mcMap);

    /**
     *  @brief  Visualise the MC particles according to their creation process
     *
     *  @param  mcMap The map from MC particles to calo hits
     **/
    void VisualiseByProcess(const LArMCParticleHelper::MCContributionMap &mcMap);

    /**
     *  @brief  Selects the MC particles to consider based on reconstructability criteria
     *
     *  @param  pMCList The input MC particle list
     *  @param  calotHitList The input calo hit list
     *  @param  mcMap The output map from MC particles to calo hits
     **/
    void MakeSelection(const pandora::MCParticleList *pMCList, const pandora::CaloHitList *pCaloHitList,
        LArMCParticleHelper::MCContributionMap &mcMap);

    std::string     m_caloHitListName;    ///< Name of input calo hit list
    bool            m_foldToPrimaries;    ///< Whether to fold the hierarchy back to primaries
    bool            m_visualisePdg;       ///< Whether or not to visualise according to the particle id
    bool            m_visualiseProcess;   ///< Whether or not to visualise according to the process that created the particle
};

} // namespace lar_content

#endif // LAR_MC_VISUAL_MONITORING_ALGORITHM_H

