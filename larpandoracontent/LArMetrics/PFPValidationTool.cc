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
    m_pNuVertexList(nullptr),
    m_nuVertexListName("NeutrinoVertices3D"),
    m_eventNumber(-1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu, const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, 
    const MCParticleVector &targetMC, const PfoVector &bestRecoMatch)
{
    ++m_eventNumber;

    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Reco neutrino vertex (we'll handle if not found)
    PandoraContentApi::GetList(*pAlgorithm, m_nuVertexListName, m_pNuVertexList);

    PFPTreeVars pfpTreeVars;
    pfpTreeVars.m_run = this->GetPandora().GetRun();
    pfpTreeVars.m_subrun = this->GetPandora().GetSubrun();
    pfpTreeVars.m_event = this->GetPandora().GetEvent();

    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        const MCParticle *const pMC(targetMC.at(i));
        const Pfo *const pBestMatch(bestRecoMatch.at(i));

        this->GetMCParticleInfo(pMC, pfpTreeVars);
        this->GetMatchingInfo(mcMatchesVec, pMC, pBestMatch, pfpTreeVars);
        this->LengthValidation(pAlgorithm, pMCNu, pMC, pBestMatch, pfpTreeVars);
        this->PIDValidation(pAlgorithm, pMC, pBestMatch, pfpTreeVars);
    }

    this->FillTree(pfpTreeVars);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::GetMCParticleInfo(const MCParticle *const pMCTarget, PFPTreeVars &pfpTreeVars)
{
    pfpTreeVars.m_trueEnergy.push_back(pMCTarget->GetEnergy());
    const CartesianVector &trueEnd(pMCTarget->GetEndpoint());
    pfpTreeVars.m_trueEndX.push_back(trueEnd.GetX());
    pfpTreeVars.m_trueEndY.push_back(trueEnd.GetY());
    pfpTreeVars.m_trueEndZ.push_back(trueEnd.GetZ());

    // Angles
    const CartesianVector &mcMom(pMCTarget->GetMomentum());

    if (mcMom.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
    {
        pfpTreeVars.m_trueThetaXZ.push_back(-4.f);
        pfpTreeVars.m_trueThetaYZ.push_back(-4.f);
        pfpTreeVars.m_trueDirX.push_back(-9999.f);
        pfpTreeVars.m_trueDirY.push_back(-9999.f);
        pfpTreeVars.m_trueDirZ.push_back(-9999.f);
    }
    else
    {
        pfpTreeVars.m_trueThetaXZ.push_back(atan2(mcMom.GetX(), mcMom.GetZ()));
        pfpTreeVars.m_trueThetaYZ.push_back(asin(mcMom.GetY() / mcMom.GetMagnitude()));
        const CartesianVector trueDir(mcMom.GetUnitVector());
        pfpTreeVars.m_trueDirX.push_back(trueDir.GetX());
        pfpTreeVars.m_trueDirY.push_back(trueDir.GetY());
        pfpTreeVars.m_trueDirZ.push_back(trueDir.GetZ());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::GetMatchingInfo(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticle *const pMCTarget, 
    const Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars)
{
    // Sorry for looping over matches again :( 
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() != pMCTarget)
            continue;

        const CaloHitList &mcHits(mcMatches.GetMC()->GetCaloHits());
        pfpTreeVars.m_nMCHits2D.push_back(mcHits.size());
        
        for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            // Get vectors to fill
            IntVector &viewNMCHits(hitType == TPC_VIEW_U ? pfpTreeVars.m_nMCHitsU : hitType == TPC_VIEW_V ? pfpTreeVars.m_nMCHitsV : 
                pfpTreeVars.m_nMCHitsW);
            IntVector &viewNPfoHits(hitType == TPC_VIEW_U ? pfpTreeVars.m_nPfoHitsU : hitType == TPC_VIEW_V ? pfpTreeVars.m_nPfoHitsV : 
                pfpTreeVars.m_nPfoHitsW);
            FloatVector &viewCompleteness(hitType == TPC_VIEW_U ? pfpTreeVars.m_completenessU : hitType == TPC_VIEW_V ? pfpTreeVars.m_completenessV : 
                pfpTreeVars.m_completenessW);
            FloatVector &viewPurity(hitType == TPC_VIEW_U ? pfpTreeVars.m_purityU : hitType == TPC_VIEW_V ? pfpTreeVars.m_purityV : 
                pfpTreeVars.m_purityW);

            // Calculate!
            float totalEnergy(0.f);
            CaloHitVector viewMCHits;
            this->GetHitsOfType(mcHits, hitType, viewMCHits, totalEnergy);
            
            viewNMCHits.push_back(viewMCHits.size());
            
            if (pBestMatch)
            {
                CaloHitList viewPfoHits;
                LArPfoHelper::GetCaloHits(pBestMatch, hitType, viewPfoHits);
                viewNPfoHits.push_back(viewPfoHits.size());
                const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode(mcMatches.GetRecoMatches().front());
                viewCompleteness.push_back(mcMatches.GetCompleteness(pRecoNode, hitType));
                viewPurity.push_back(mcMatches.GetPurity(pRecoNode, hitType));
            }
            else
            {
                viewNPfoHits.push_back(0);
                viewCompleteness.push_back(0);
                viewPurity.push_back(0);
            }
        }
             
        if (pBestMatch)
        {
            const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode(mcMatches.GetRecoMatches().front());
            pfpTreeVars.m_hasMatch.push_back(1);
            pfpTreeVars.m_nPfoHits2D.push_back(LArPfoHelper::GetNumberOfTwoDHits(pBestMatch));
            pfpTreeVars.m_nPfoHits3D.push_back(LArPfoHelper::GetNumberOfThreeDHits(pBestMatch));
            pfpTreeVars.m_completeness.push_back(mcMatches.GetCompleteness(pRecoNode));
            pfpTreeVars.m_purity.push_back(mcMatches.GetPurity(pRecoNode));
        }
        else
        {
            pfpTreeVars.m_hasMatch.push_back(0);
            pfpTreeVars.m_nPfoHits2D.push_back(-1);
            pfpTreeVars.m_nPfoHits3D.push_back(-1);
            pfpTreeVars.m_completeness.push_back(0.f);
            pfpTreeVars.m_purity.push_back(0.f);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::LengthValidation(const Algorithm *const pAlgorithm, const MCParticle *const pMCNu, 
    const MCParticle *const pMCTarget, const Pfo *const pBestMatch, PFPTreeVars &pfpTreeVars)
{
    // True neutrino vertex
    const CartesianVector &trueNuVertex(pMCNu->GetVertex());

    // Fill
    const CartesianVector &trueVertex(pMCTarget->GetVertex());
    pfpTreeVars.m_trueVertexX.push_back(trueVertex.GetX());
    pfpTreeVars.m_trueVertexY.push_back(trueVertex.GetY());
    pfpTreeVars.m_trueVertexZ.push_back(trueVertex.GetZ());
    pfpTreeVars.m_trueLength.push_back((trueVertex - pMCTarget->GetEndpoint()).GetMagnitude());
    pfpTreeVars.m_trueDisplacement.push_back((trueNuVertex - trueVertex).GetMagnitude());

    if (!pBestMatch)
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
        const Vertex *pRecoVertex(nullptr);

        try
        {
            pRecoVertex = LArPfoHelper::GetVertex(pBestMatch);
        }
        catch(...){}

        if (pRecoVertex)
        {
            pfpTreeVars.m_recoVertexX.push_back(pRecoVertex->GetPosition().GetX());
            pfpTreeVars.m_recoVertexY.push_back(pRecoVertex->GetPosition().GetY());
            pfpTreeVars.m_recoVertexZ.push_back(pRecoVertex->GetPosition().GetZ());

            // Signed vertexAcc
            const float vertexAcc((pRecoVertex->GetPosition() - trueVertex).GetMagnitude());
            const float sign((vertexAcc < std::numeric_limits<float>::epsilon() || pMCTarget->GetEnergy() < std::numeric_limits<float>::epsilon()) ? 1.f : 
                             (pRecoVertex->GetPosition() - trueVertex).GetOpeningAngle(pMCTarget->GetMomentum()) < (M_PI * 0.5) ? 1.f : -1.f);
            pfpTreeVars.m_vertexAcc.push_back(vertexAcc * sign);
            
            try
            {
                pfpTreeVars.m_recoLength.push_back(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pBestMatch))); // TODO - better fit.
            }
            catch(...)
            {
                pfpTreeVars.m_recoLength.push_back(-1.f);
            }

            if (m_pNuVertexList && !m_pNuVertexList->empty())
            {
                pfpTreeVars.m_recoDisplacement.push_back((m_pNuVertexList->front()->GetPosition() - pRecoVertex->GetPosition()).GetMagnitude());
            }
            else
            {
                pfpTreeVars.m_recoDisplacement.push_back(-1.f);
            }
        }
        else
        {
            pfpTreeVars.m_recoVertexX.push_back(-9999.f);
            pfpTreeVars.m_recoVertexY.push_back(-9999.f);
            pfpTreeVars.m_recoVertexZ.push_back(-9999.f);
            pfpTreeVars.m_vertexAcc.push_back(-9999.f);
            pfpTreeVars.m_recoLength.push_back(-1.f);
            pfpTreeVars.m_recoDisplacement.push_back(-1.f);
        }    
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::PIDValidation(const Algorithm *const pAlgorithm, const MCParticle *const pMCTarget, const Pfo *const pBestMatch, 
    PFPTreeVars &pfpTreeVars)
{
    pfpTreeVars.m_truePDG.push_back(pMCTarget->GetParticleId());

    if (pBestMatch)
    {
        pfpTreeVars.m_isTrack.push_back(LArPfoHelper::IsTrack(pBestMatch));
        pfpTreeVars.m_isShower.push_back(LArPfoHelper::IsShower(pBestMatch));
    }
    else
    {
        pfpTreeVars.m_isTrack.push_back(-1);
        pfpTreeVars.m_isShower.push_back(-1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFPValidationTool::FillTree(PFPTreeVars &pfpTreeVars)
{
    IntVector &truePDG = pfpTreeVars.m_truePDG;
    FloatVector &trueEnergy = pfpTreeVars.m_trueEnergy;
    FloatVector &trueThetaXZ = pfpTreeVars.m_trueThetaXZ;
    FloatVector &trueThetaYZ = pfpTreeVars.m_trueThetaYZ;
    IntVector &isTrack = pfpTreeVars.m_isTrack;
    IntVector &isShower = pfpTreeVars.m_isShower;
    IntVector &hasMatch = pfpTreeVars.m_hasMatch;
    IntVector &nMCHitsU = pfpTreeVars.m_nMCHitsU;
    IntVector &nMCHitsV = pfpTreeVars.m_nMCHitsV;
    IntVector &nMCHitsW = pfpTreeVars.m_nMCHitsW;
    IntVector &nMCHits2D = pfpTreeVars.m_nMCHits2D;
    IntVector &nPfoHitsU = pfpTreeVars.m_nPfoHitsU;
    IntVector &nPfoHitsV = pfpTreeVars.m_nPfoHitsV;
    IntVector &nPfoHitsW = pfpTreeVars.m_nPfoHitsW;
    IntVector &nPfoHits2D = pfpTreeVars.m_nPfoHits2D;
    IntVector &nPfoHits3D = pfpTreeVars.m_nPfoHits3D;
    FloatVector &completeness = pfpTreeVars.m_completeness;
    FloatVector &completenessU = pfpTreeVars.m_completenessU;
    FloatVector &completenessV = pfpTreeVars.m_completenessV;
    FloatVector &completenessW = pfpTreeVars.m_completenessW;
    FloatVector &purity = pfpTreeVars.m_purity;
    FloatVector &purityU = pfpTreeVars.m_purityU;
    FloatVector &purityV = pfpTreeVars.m_purityV;
    FloatVector &purityW = pfpTreeVars.m_purityW;
    FloatVector &trueVertexX = pfpTreeVars.m_trueVertexX;
    FloatVector &trueVertexY = pfpTreeVars.m_trueVertexY; 
    FloatVector &trueVertexZ = pfpTreeVars.m_trueVertexZ;
    FloatVector &trueEndX = pfpTreeVars.m_trueEndX;
    FloatVector &trueEndY = pfpTreeVars.m_trueEndY; 
    FloatVector &trueEndZ = pfpTreeVars.m_trueEndZ;
    FloatVector &trueDirX = pfpTreeVars.m_trueDirX;
    FloatVector &trueDirY = pfpTreeVars.m_trueDirY; 
    FloatVector &trueDirZ = pfpTreeVars.m_trueDirZ;
    FloatVector &trueLength = pfpTreeVars.m_trueLength;
    FloatVector &trueDisplacement = pfpTreeVars.m_trueDisplacement;
    FloatVector &recoVertexX = pfpTreeVars.m_recoVertexX;
    FloatVector &recoVertexY = pfpTreeVars.m_recoVertexY;
    FloatVector &recoVertexZ = pfpTreeVars.m_recoVertexZ;
    FloatVector &vertexAcc = pfpTreeVars.m_vertexAcc;
    FloatVector &recoLength = pfpTreeVars.m_recoLength;
    FloatVector &recoDisplacement = pfpTreeVars.m_recoDisplacement;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "Run", pfpTreeVars.m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "Subrun", pfpTreeVars.m_subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "Event", pfpTreeVars.m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "EventCount", m_eventNumber));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TruePDG", &truePDG));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TrueEnergy", &trueEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TrueThetaXZ", &trueThetaXZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_TrueThetaYZ", &trueThetaYZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_IsTrack", &isTrack));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_IsShower", &isShower));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_HasMatch", &hasMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHitsU", &nMCHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHitsV", &nMCHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHitsW", &nMCHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_NMCHits2D", &nMCHits2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHitsU", &nPfoHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHitsV", &nPfoHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHitsW", &nPfoHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHits2D", &nPfoHits2D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_NPfoHits3D", &nPfoHits3D));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Completeness", &completeness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessU", &completenessU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessV", &completenessV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_CompletenessW", &completenessW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Purity", &purity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityU", &purityU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityV", &purityV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_PurityW", &purityW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_VertexX", &trueVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_VertexY", &trueVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_VertexZ", &trueVertexZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndX", &trueEndX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndY", &trueEndY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_EndZ", &trueEndZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_DirX", &trueDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_DirY", &trueDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_DirZ", &trueDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_Length", &trueLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "MCP_Displacement", &trueDisplacement));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexX", &recoVertexX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexY", &recoVertexY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexZ", &recoVertexZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_VertexAcc", &vertexAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Length", &recoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFPTree", "BM_Displacement", &recoDisplacement));
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
