/**
 *  @file   larpandoracontent/LArMetrics/EventValidationTool.h
 *
 *  @brief  Header file for the track validation tool class.
 *
 *  $Log: $
 */
#ifndef EVENT_VALIDATION_TOOL_H
#define EVENT_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  EventValidationTool class
 */
class EventValidationTool : public BaseValidationTool
{
public:

struct EventTreeVars
{

    EventTreeVars();

    int m_run;
    int m_subrun;
    int m_event;
    int m_nTargets;
    int m_nuPDG;
    float m_nuEnergy;
    float m_nuVisEnergy;
    int m_isCC;
    pandora::CartesianVector m_trueNuVertex;
    pandora::CartesianVector m_recoNuVertexPass1;
    pandora::CartesianVector m_recoNuVertexPass2;
    float m_recoNuVertexAccPass1;
    float m_recoNuVertexAccPass2;
};
    /**
     *  @brief  Default constructor
     */
    EventValidationTool();

    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetNeutrinoVariables(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, EventTreeVars &eventTreeVars);
    void GetInteractionTypeVariables(const pandora::MCParticle *const pMCNu, EventTreeVars &eventTreeVars);
    void FillTree(EventTreeVars &eventTreeVars);

    std::string m_nuVertexPass1ListName;
    std::string m_nuVertexPass2ListName;
    int m_eventNumber;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline EventValidationTool::EventTreeVars::EventTreeVars() :
    m_run(-1),
    m_subrun(-1),
    m_event(-1),
    m_nTargets(-1),
    m_nuPDG(-1),
    m_nuEnergy(-1.f),
    m_nuVisEnergy(-1.f),
    m_isCC(-1),
    m_trueNuVertex(pandora::CartesianVector(-9999.f, -9999.f, -9999.f)),
    m_recoNuVertexPass1(pandora::CartesianVector(-9999.f, -9999.f, -9999.f)),
    m_recoNuVertexPass2(pandora::CartesianVector(-9999.f, -9999.f, -9999.f)),
    m_recoNuVertexAccPass1(-1.f),
    m_recoNuVertexAccPass2(-1.f)
{
}

} // namespace lar_content

#endif // #ifndef EVENT_VALIDATION_TOOL_H
