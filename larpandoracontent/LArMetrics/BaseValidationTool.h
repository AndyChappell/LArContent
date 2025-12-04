/**
 *  @file   larpandoracontent/LArMetrics/BaseValidationTool.h
 *
 *  @brief  Header file for the pfp validation tool class.
 *
 *  $Log: $
 */
#ifndef BASE_VALIDATION_TOOL_H
#define BASE_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

namespace lar_content
{

/**
 *  @brief  BaseValidationTool class
 */
class BaseValidationTool : public pandora::AlgorithmTool
{
public:

    virtual void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch) = 0;

protected:
    void GetHitsOfType(const pandora::CaloHitList &inputList, const pandora::HitType hitType, pandora::CaloHitVector &outputVector, 
        float &totalEnergy);
};

} // namespace lar_content

#endif // #ifndef BASE_VALIDATION_TOOL_H
