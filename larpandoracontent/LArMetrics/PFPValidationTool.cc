/**
 *  @file   larpandoracontent/LArMetrics/PFPValidationTool.cc
 *
 *  @brief  Implementation of the pfp validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMetrics/PFPValidationTool.h"

using namespace pandora;

namespace lar_content
{

PFPValidationTool::PFPValidationTool() :
    m_nuVertexListName("NeutrinoVertices3D")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu, const MCParticleVector &targetMC, 
    const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    PFPTreeVars pfpTreeVars;
    this->LengthValidation(pAlgorithm, pMCNu, targetMC, bestRecoMatch, pfpTreeVars);
    this->PIDValidation(pAlgorithm, targetMC, bestRecoMatch, pfpTreeVars);
    this->FillTree(pfpTreeVars);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::LengthValidation(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu, 
    const MCParticleVector &targetMC, const PfoVector &bestRecoMatch, PFPTreeVars &pfpTreeVars)
{
    // True neutrino vertex
    const CartesianVector &trueNuVertex(pMCNu->GetVertex());

    // Reco neutrino vertex
    const VertexList *pNuVertexList(nullptr);
    PandoraContentApi::GetList(*pAlgorithm, m_nuVertexListName, pNuVertexList);

    // Fill
    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        const MCParticle *const pMCParticle(targetMC.at(i));
        const CartesianVector &trueVertex(pMCParticle->GetVertex());
        pfpTreeVars.m_trueVertexX.push_back(trueVertex.GetX());
        pfpTreeVars.m_trueVertexY.push_back(trueVertex.GetY());
        pfpTreeVars.m_trueVertexZ.push_back(trueVertex.GetZ());
        pfpTreeVars.m_trueLength.push_back((trueVertex - pMCParticle->GetEndpoint()).GetMagnitude());
        pfpTreeVars.m_trueDisplacement.push_back((trueNuVertex - trueVertex).GetMagnitude());

        if (!bestRecoMatch.at(i))
        {
            pfpTreeVars.m_recoVertexX.push_back(-9999.f);
            pfpTreeVars.m_recoVertexY.push_back(-9999.f);
            pfpTreeVars.m_recoVertexZ.push_back(-9999.f);
            pfpTreeVars.m_vertexAcc.push_back(-9999.f);
            pfpTreeVars.m_recoLength.push_back(-1.f);
            pfpTreeVars.m_recoDisplacement.push_back(-1.f);
        }
        else
        {
            const Pfo *const pPfo(bestRecoMatch.at(i));
            const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(pPfo));

            pfpTreeVars.m_recoVertexX.push_back(pRecoVertex->GetPosition().GetX());
            pfpTreeVars.m_recoVertexY.push_back(pRecoVertex->GetPosition().GetY());
            pfpTreeVars.m_recoVertexZ.push_back(pRecoVertex->GetPosition().GetZ());
            pfpTreeVars.m_vertexAcc.push_back((pRecoVertex->GetPosition() - trueVertex).GetMagnitude()); // TODO - sign this // TODO - Add wrong end?

            try
            {
                pfpTreeVars.m_recoLength.push_back(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pPfo))); // TODO - better fit.
            }
            catch(...)
            {
                pfpTreeVars.m_recoLength.push_back(-1.f);
            }

            if (pNuVertexList && !pNuVertexList->empty())
            {
                pfpTreeVars.m_recoDisplacement.push_back((pNuVertexList->front()->GetPosition() - pRecoVertex->GetPosition()).GetMagnitude());
            }
            else
            {
                pfpTreeVars.m_recoDisplacement.push_back(-1.f);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::PIDValidation(const Algorithm *const pAlgorithm, const MCParticleVector &targetMC, const PfoVector &bestRecoMatch, 
    PFPTreeVars &pfpTreeVars)
{
    for (unsigned int iMC = 0; iMC < targetMC.size(); ++iMC)
    {
        const MCParticle *const pMCTarget(targetMC.at(iMC));

        pfpTreeVars.m_truePDG.push_back(pMCTarget->GetParticleId());

        if (bestRecoMatch.at(iMC))
        {
            const Pfo *const pRecoMatch(bestRecoMatch.at(iMC));
            pfpTreeVars.m_isTrack.push_back(LArPfoHelper::IsTrack(pRecoMatch));
            pfpTreeVars.m_isShower.push_back(LArPfoHelper::IsShower(pRecoMatch));
        }
        else
        {
            pfpTreeVars.m_isTrack.push_back(-1);
            pfpTreeVars.m_isShower.push_back(-1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::FillTree(PFPTreeVars &pfpTreeVars)
{
    IntVector& truePDG = pfpTreeVars.m_truePDG;
    IntVector& isTrack = pfpTreeVars.m_isTrack;
    IntVector& isShower = pfpTreeVars.m_isShower;
    FloatVector& trueVertexX = pfpTreeVars.m_trueVertexX;
    FloatVector& trueVertexY = pfpTreeVars.m_trueVertexY; 
    FloatVector& trueVertexZ = pfpTreeVars.m_trueVertexZ;
    FloatVector& trueLength = pfpTreeVars.m_trueLength;
    FloatVector& trueDisplacement = pfpTreeVars.m_trueDisplacement;
    FloatVector& recoVertexX = pfpTreeVars.m_recoVertexX;
    FloatVector& recoVertexY = pfpTreeVars.m_recoVertexY;
    FloatVector& recoVertexZ = pfpTreeVars.m_recoVertexZ;
    FloatVector& vertexAcc = pfpTreeVars.m_vertexAcc;
    FloatVector& recoLength = pfpTreeVars.m_recoLength;
    FloatVector& recoDisplacement = pfpTreeVars.m_recoDisplacement;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TruePDG", &truePDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_IsTrack", &isTrack));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_IsShower", &isShower));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_VertexX", &trueVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_VertexY", &trueVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_VertexZ", &trueVertexZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_Length", &trueLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_Displacement", &trueDisplacement));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_VertexX", &recoVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_VertexY", &recoVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_VertexZ", &recoVertexZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_VertexAcc", &vertexAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_Length", &recoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_Displacement", &recoDisplacement));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "PFPTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PFPValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexListName", m_nuVertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
