/**
 *  @file   larpandoracontent/LArMetrics/ShowerValidationTool.h
 *
 *  @brief  Header file for the shower validation tool class.
 *
 *  $Log: $
 */
#ifndef SHOWER_VALIDATION_TOOL_H
#define SHOWER_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  ShowerValidationTool class
 */
class ShowerValidationTool : public BaseValidationTool
{
public:

struct ShowerTreeVars
{
    int m_run;
    int m_subrun;
    int m_event;
    pandora::FloatVector m_trueShrDirX;
    pandora::FloatVector m_trueShrDirY;
    pandora::FloatVector m_trueShrDirZ;
    pandora::FloatVector m_coreTrueLengthFromU;
    pandora::FloatVector m_coreTrueLengthFromV;
    pandora::FloatVector m_coreTrueLengthFromW;
    pandora::FloatVector m_recoShrVtxX;
    pandora::FloatVector m_recoShrVtxY;
    pandora::FloatVector m_recoShrVtxZ;
    pandora::FloatVector m_recoShrDirX;
    pandora::FloatVector m_recoShrDirY;
    pandora::FloatVector m_recoShrDirZ;
    pandora::FloatVector m_coreRecoLength;
    pandora::FloatVector m_recoShrLength;
    pandora::FloatVector m_recoShrDirAcc;
    pandora::FloatVector m_moliereRadius;
    pandora::IntVector m_nInitialMCHits;
    pandora::IntVector m_nInitialMCHitsU;
    pandora::IntVector m_nInitialMCHitsV;
    pandora::IntVector m_nInitialMCHitsW;
    pandora::IntVector m_nInitialPfoHits;
    pandora::IntVector m_nInitialPfoHitsU;
    pandora::IntVector m_nInitialPfoHitsV;
    pandora::IntVector m_nInitialPfoHitsW;
    pandora::FloatVector m_initialCompleteness;
    pandora::FloatVector m_initialCompletenessU;
    pandora::FloatVector m_initialCompletenessV;
    pandora::FloatVector m_initialCompletenessW;
    pandora::FloatVector m_initialPurity;
    pandora::FloatVector m_initialPurityU;
    pandora::FloatVector m_initialPurityV;
    pandora::FloatVector m_initialPurityW;
};

    /**
     *  @brief  Default constructor
     */
    ShowerValidationTool();

    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

    bool FitShower(const pandora::Pfo *const pPfo, pandora::CartesianVector &showerVertex, 
        pandora::CartesianVector &showerDirection, float &showerLength);

    void GetMoliere(const pandora::Pfo *const pPfo, const pandora::CartesianVector &showerVertex, const pandora::CartesianVector &showerDirection,
        ShowerTreeVars &showerTreeVars);

    void GetTrueLength(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticle *const pMCParticle,
        ShowerTreeVars &showerTreeVars);

    void GetHitsOfType(const pandora::CaloHitList &inputList, const pandora::HitType hitType, pandora::CaloHitVector &outputVector, float &totalEnergy);

    void GetInitialRegionVars(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticle *const pMCParticle, 
        const pandora::Pfo *const pPfo, ShowerTreeVars &showerTreeVars);

    void FillTree(ShowerTreeVars &showerTreeVars);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_trueLengthEnergyFrac;
    float m_initialRegion3D;
};

} // namespace lar_content

#endif // #ifndef SHOWER_VALIDATION_TOOL_H
