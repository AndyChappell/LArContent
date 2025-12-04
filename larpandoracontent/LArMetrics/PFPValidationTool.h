/**
 *  @file   larpandoracontent/LArMetrics/PFPValidationTool.h
 *
 *  @brief  Header file for the pfp validation tool class.
 *
 *  $Log: $
 */
#ifndef PFP_VALIDATION_TOOL_H
#define PFP_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  PFPValidationTool class
 */
class PFPValidationTool : public BaseValidationTool
{
public:

struct PFPTreeVars
{
    int m_run;
    int m_subrun;
    int m_event;
    pandora::IntVector m_truePDG;
    pandora::IntVector m_isTrack;
    pandora::IntVector m_isShower;
    pandora::IntVector m_nMCHitsU;
    pandora::IntVector m_nMCHitsV;
    pandora::IntVector m_nMCHitsW;
    pandora::IntVector m_nPfoHitsU;
    pandora::IntVector m_nPfoHitsV;
    pandora::IntVector m_nPfoHitsW;
    pandora::FloatVector m_completeness;
    pandora::FloatVector m_completenessU;
    pandora::FloatVector m_completenessV;
    pandora::FloatVector m_completenessW;
    pandora::FloatVector m_purity;
    pandora::FloatVector m_purityU;
    pandora::FloatVector m_purityV;
    pandora::FloatVector m_purityW;
    pandora::FloatVector m_trueVertexX;
    pandora::FloatVector m_trueVertexY; 
    pandora::FloatVector m_trueVertexZ;
    pandora::FloatVector m_trueLength;
    pandora::FloatVector m_trueDisplacement;
    pandora::FloatVector m_recoVertexX;
    pandora::FloatVector m_recoVertexY;
    pandora::FloatVector m_recoVertexZ;
    pandora::FloatVector m_vertexAcc;
    pandora::FloatVector m_recoLength;
    pandora::FloatVector m_recoDisplacement;
};

    /**
     *  @brief  Default constructor
     */
    PFPValidationTool();

    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

    void GetMatchingInfo(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticle *const pMCTarget, 
        const pandora::Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars);

    void LengthValidation(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const pandora::MCParticle *const pTargetMC, const pandora::Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars);

    void PIDValidation(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pTargetMC, 
        const pandora::Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars);

    void FillTree(PFPTreeVars &pfpTreeVars);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    const pandora::VertexList *m_pNuVertexList;

    std::string m_nuVertexListName;
    int m_eventNumber;
};

} // namespace lar_content

#endif // #ifndef PFP_VALIDATION_TOOL_H
