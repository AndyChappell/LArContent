/**
 *  @file   larpandoracontent/LArMetrics/HierarchyValidationTool.h
 *
 *  @brief  Header file for the hierarchy validation tool class.
 *
 *  $Log: $
 */
#ifndef HIERARCHY_VALIDATION_TOOL_H
#define HIERARCHY_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  HierarchyValidationTool class
 */
class HierarchyValidationTool : public BaseValidationTool
{

struct HierarchyTreeVars
{
    int m_run;
    int m_subrun;
    int m_event;
    pandora::IntVector m_trueTier;
    pandora::IntVector m_trueParentIndex;
    pandora::IntVector m_recoTier;
    pandora::IntVector m_recoParentIndex;
};

public:
    /**
     *  @brief  Default constructor
     */
    HierarchyValidationTool();

    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

private:
    typedef std::map<const pandora::MCParticle*, std::pair<const pandora::MCParticle*, int>> Hierarchy;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void BuildVisibleHierarchy(const pandora::MCParticle *const pMCParticle, const pandora::MCParticle *const pMCParent, 
        const pandora::MCParticleVector &targetMC, const int childTier, Hierarchy &hierarchy);

    void FillTrueVariables(const pandora::MCParticle *const pMC, const pandora::MCParticleVector &targetMC,
        const Hierarchy &hierarchy, HierarchyTreeVars &hierarchyTreeVars);

    void FillRecoVariables(const pandora::Pfo *const pBestMatch, const pandora::PfoVector &bestMatches, 
        HierarchyTreeVars &hierarchyTreeVars);

    void FillTree(HierarchyTreeVars &hierarchyTreeVars);
};

} // namespace lar_content

#endif // #ifndef HIERARCHY_VALIDATION_TOOL_H
