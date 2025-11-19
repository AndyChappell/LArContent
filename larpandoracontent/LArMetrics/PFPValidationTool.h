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

namespace lar_content
{

/**
 *  @brief  PFPValidationTool class
 */
class PFPValidationTool : public pandora::AlgorithmTool
{
public:

struct PFPTreeVars
{
    pandora::IntVector m_truePDG;
    pandora::IntVector m_isTrack;
    pandora::IntVector m_isShower;
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
        const pandora::MCParticleVector &targetMC, const pandora::PfoVector &bestRecoMatch);

    void LengthValidation(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const pandora::MCParticleVector &targetMC, const pandora::PfoVector &bestRecoMatch, PFPTreeVars &pfpTreeVars);

    void PIDValidation(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch, PFPTreeVars &pfpTreeVars);

    void FillTree(PFPTreeVars &pfpTreeVars);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_nuVertexListName;
};

} // namespace lar_content

#endif // #ifndef PFP_VALIDATION_TOOL_H
