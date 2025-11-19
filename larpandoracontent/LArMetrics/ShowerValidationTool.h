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

namespace lar_content
{

/**
 *  @brief  ShowerValidationTool class
 */
class ShowerValidationTool : public pandora::AlgorithmTool
{
public:

struct ShowerTreeVars
{
    pandora::FloatVector m_trueShrDirX;
    pandora::FloatVector m_trueShrDirY;
    pandora::FloatVector m_trueShrDirZ;
    pandora::FloatVector m_recoShrVtxX;
    pandora::FloatVector m_recoShrVtxY;
    pandora::FloatVector m_recoShrVtxZ;
    pandora::FloatVector m_recoShrDirX;
    pandora::FloatVector m_recoShrDirY;
    pandora::FloatVector m_recoShrDirZ;
    pandora::FloatVector m_recoShrDirAcc;
    pandora::FloatVector m_moliereRadius;
};

    /**
     *  @brief  Default constructor
     */
    ShowerValidationTool();

    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

    void GetShowerDirection(const pandora::Pfo *const pPfo, pandora::CartesianVector &showerVertex, 
        pandora::CartesianVector &showerDirection);

    void GetMoliere(const pandora::Pfo *const pPfo, const pandora::CartesianVector &showerVertex, const pandora::CartesianVector &showerDirection,
        ShowerTreeVars &showerTreeVars);

    void FillTree(ShowerTreeVars &showerTreeVars);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef SHOWER_VALIDATION_TOOL_H
