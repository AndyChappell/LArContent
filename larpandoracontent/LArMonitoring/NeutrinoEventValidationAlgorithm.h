/**
 *  @file   larpandoracontent/LArMonitoring/NeutrinoEventValidationAlgorithm.h
 *
 *  @brief  Header file for the neutrino event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_EVENT_VALIDATION_ALGORITHM_H
#define LAR_NEUTRINO_EVENT_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationBaseAlgorithm.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

#include <map>

namespace lar_content
{

/**
 *  @brief  NeutrinoEventValidationAlgorithm class
 */
class NeutrinoEventValidationAlgorithm : public EventValidationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoEventValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~NeutrinoEventValidationAlgorithm();

private:
    /**
     *  @brief  Fill the validation info containers
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  validationInfo to receive the validation info
     */
    void FillValidationInfo(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList,
        const pandora::PfoList *const pPfoList, ValidationInfo &validationInfo) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, unsigned int> PfoToIdMap;

    /**
     *  @brief  Print matching information in a provided validation info object, and write information to tree if configured to do so
     *
     *  @param  validationInfo the validation info
     *  @param  useInterpretedMatching whether to use the interpreted (rather than raw) matching information
     *  @param  printToScreen whether to print the information to screen
     *  @param  fillTree whether to write the information to tree
     */
    void ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<pandora::HitType> HitTypeVector;

    bool m_useTrueNeutrinosOnly; ///< Whether to consider only mc particles that were neutrino induced
};

} // namespace lar_content

#endif // LAR_NEUTRINO_EVENT_VALIDATION_ALGORITHM_H
