/**
 *  @file   larpandoracontent/LArMetrics/TrackValidationTool.h
 *
 *  @brief  Header file for the track validation tool class.
 *
 *  $Log: $
 */
#ifndef TRACK_VALIDATION_TOOL_H
#define TRACK_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  TrackValidationTool class
 */
class TrackValidationTool : public BaseValidationTool
{

struct TrackTreeVars
{
    int m_run;
    int m_subrun;
    int m_event;
    pandora::IntVector m_hasMichel; 
    pandora::IntVector m_michelFromMuon; 
    pandora::IntVector m_hasTargetMichel;
    pandora::IntVector m_hasRecoMichel;
    pandora::IntVector m_michelIndex; 
    pandora::IntVector m_michelIsChild;
    pandora::IntVector m_michelIsShower;
};

public:
    /**
     *  @brief  Default constructor
     */
    TrackValidationTool();

    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

    void MichelValidation(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch, TrackTreeVars &trackTreeVars);

    void FillTree(TrackTreeVars &trackTreeVars);
private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef TRACK_VALIDATION_TOOL_H
