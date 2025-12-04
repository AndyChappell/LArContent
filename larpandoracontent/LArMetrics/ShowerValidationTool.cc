/**
 *  @file   larpandoracontent/LArMetrics/ShowerValidationTool.cc
 *
 *  @brief  Implementation of the shower validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/ShowerValidationTool.h"

using namespace pandora;

namespace lar_content
{

ShowerValidationTool::ShowerValidationTool() :
    m_trueLengthEnergyFrac(0.95),
    m_initialRegion3D(14.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const /*pMCNu*/, 
    const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticleVector &targetMC, 
    const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ShowerTreeVars showerTreeVars;
    showerTreeVars.m_run = this->GetPandora().GetRun();
    showerTreeVars.m_subrun = this->GetPandora().GetSubrun();
    showerTreeVars.m_event = this->GetPandora().GetEvent();

    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        // Some truth variables
        const MCParticle *const pMC(targetMC.at(i));
        const CartesianVector trueDir(pMC->GetMomentum().GetUnitVector());
        showerTreeVars.m_trueShrDirX.push_back(trueDir.GetX());
        showerTreeVars.m_trueShrDirY.push_back(trueDir.GetY());
        showerTreeVars.m_trueShrDirZ.push_back(trueDir.GetZ());
        this->GetTrueLength(mcMatchesVec, pMC, showerTreeVars);

        // Look at reco 
        const Pfo *const pPfo(bestRecoMatch.at(i));
        bool success(pPfo);

        if (success)
        {
            float recoShrLength(-1.f);
            CartesianVector recoShrVtx(0.f, 0.f, 0.f), recoShrDir(0.f, 0.f, 0.f);
            success = this->FitShower(pPfo, recoShrVtx, recoShrDir, recoShrLength);

            if (success)
            {
                showerTreeVars.m_recoShrVtxX.push_back(recoShrVtx.GetX());
                showerTreeVars.m_recoShrVtxY.push_back(recoShrVtx.GetY());
                showerTreeVars.m_recoShrVtxZ.push_back(recoShrVtx.GetZ());
                showerTreeVars.m_recoShrDirX.push_back(recoShrDir.GetX());
                showerTreeVars.m_recoShrDirY.push_back(recoShrDir.GetY());
                showerTreeVars.m_recoShrDirZ.push_back(recoShrDir.GetZ());
                showerTreeVars.m_recoShrLength.push_back(recoShrLength);
                showerTreeVars.m_recoShrDirAcc.push_back(trueDir.GetOpeningAngle(recoShrDir));

                this->GetMoliere(pPfo, recoShrVtx, recoShrDir, showerTreeVars);
                this->GetInitialRegionVars(mcMatchesVec, pMC, pPfo, showerTreeVars);
            }
        }

        if (!success)
        {
            showerTreeVars.m_recoShrVtxX.push_back(-9999.f);
            showerTreeVars.m_recoShrVtxY.push_back(-9999.f);
            showerTreeVars.m_recoShrVtxZ.push_back(-9999.f);
            showerTreeVars.m_recoShrDirX.push_back(-9999.f);
            showerTreeVars.m_recoShrDirY.push_back(-9999.f);
            showerTreeVars.m_recoShrDirZ.push_back(-9999.f);
            showerTreeVars.m_recoShrLength.push_back(-1.f);
            showerTreeVars.m_recoShrDirAcc.push_back(-1.f);
            showerTreeVars.m_moliereRadius.push_back(-1.f);
            showerTreeVars.m_coreRecoLength.push_back(-1.f);
            showerTreeVars.m_initialCompletenessU.push_back(-1.f);
            showerTreeVars.m_initialCompletenessV.push_back(-1.f);
            showerTreeVars.m_initialCompletenessW.push_back(-1.f);
            showerTreeVars.m_initialPurityU.push_back(-1.f);
            showerTreeVars.m_initialPurityV.push_back(-1.f);
            showerTreeVars.m_initialPurityW.push_back(-1.f);
        }
    }

    this->FillTree(showerTreeVars);
}

//------------------------------------------------------------------------------------------------------------------------------------------

// Replicate PandoraShower fitting
bool ShowerValidationTool::FitShower(const Pfo *const pPfo, CartesianVector &showerVertex, CartesianVector &showerDirection, float &showerLength)
{
    CartesianPointVector positions3D;
    LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, positions3D);
    if (positions3D.empty()){return false;}

    const Vertex *pRecoVertex(nullptr);
    try
    {
        pRecoVertex = LArPfoHelper::GetVertex(pPfo);
    }
    catch(...){ return false;}

    const CartesianVector vertexPosition(pRecoVertex->GetPosition());
    const LArShowerPCA initialLArShowerPCA(LArPfoHelper::GetPrincipalComponents(positions3D, vertexPosition)); 

    // Ensure successful creation of all structures before placing results in output containers, remaking LArShowerPCA with updated vertex
    const CartesianVector& centroid(initialLArShowerPCA.GetCentroid());
    const CartesianVector& primaryAxis(initialLArShowerPCA.GetPrimaryAxis());
    const CartesianVector& secondaryAxis(initialLArShowerPCA.GetSecondaryAxis());
    const CartesianVector& tertiaryAxis(initialLArShowerPCA.GetTertiaryAxis());
    const CartesianVector& eigenvalues(initialLArShowerPCA.GetEigenValues());

    // Project the PFParticle vertex onto the PCA axis
    const CartesianVector projectedVertexPosition(centroid -
          primaryAxis.GetUnitVector() * (centroid - vertexPosition).GetDotProduct(primaryAxis));

    // By convention, principal axis should always point away from vertex
    const float testProjection(primaryAxis.GetDotProduct(projectedVertexPosition - centroid));
    const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

    const LArShowerPCA larShowerPCA(centroid, primaryAxis * directionScaleFactor, secondaryAxis * directionScaleFactor,
        tertiaryAxis * directionScaleFactor, eigenvalues);

    showerVertex = projectedVertexPosition;
    showerDirection = larShowerPCA.GetPrimaryAxis();
    showerLength = larShowerPCA.GetAxisLengths().GetX();

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetMoliere(const Pfo *const pPfo, const CartesianVector &showerVertex, const CartesianVector &showerDirection,
    ShowerTreeVars &showerTreeVars)
{
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHits3D);

    // Get 3D hits created from collection-view hits
    float totalEnergy(0.f);
    CaloHitVector caloHits3DFromW;

    for (const CaloHit *const pHit3D : caloHits3D)
    {
        const CaloHit *const pHit2D{static_cast<const CaloHit *>(pHit3D->GetParentAddress())};

        if (pHit2D->GetHitType() == TPC_VIEW_W)
        {
            totalEnergy += pHit2D->GetElectromagneticEnergy();
            caloHits3DFromW.push_back(pHit3D);
        }
    }

    // Now do Molliere
    float runningEnergySum(0.f), moliereRadius(-1.f);

    std::sort(caloHits3DFromW.begin(), caloHits3DFromW.end(),
        [&showerDirection, &showerVertex](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool
        {
            const CartesianVector positionA(pCaloHitA->GetPositionVector() - showerVertex);
            const CartesianVector positionB(pCaloHitB->GetPositionVector() - showerVertex);
            
            const float tA(showerDirection.GetCrossProduct(positionA).GetMagnitude());
            const float tB(showerDirection.GetCrossProduct(positionB).GetMagnitude());

            return tA < tB;
        });

    for (const CaloHit *const pHit3D : caloHits3DFromW)
    {
        const CaloHit *const pHit2D{static_cast<const CaloHit *>(pHit3D->GetParentAddress())};
        const float hitEnergy(std::fabs(pHit2D->GetElectromagneticEnergy()));
        runningEnergySum += hitEnergy;

        if ((totalEnergy > std::numeric_limits<float>::epsilon()) && 
            ((runningEnergySum / totalEnergy) > 0.9f))
        {
            const CartesianVector displacement(pHit3D->GetPositionVector() - showerVertex);
            moliereRadius = showerDirection.GetCrossProduct(displacement).GetMagnitude();
            break;
        }
    }

    showerTreeVars.m_moliereRadius.push_back(moliereRadius);

    // Now do core reco length
    runningEnergySum = 0.f;
    float coreShowerLength = 0.f;

    std::sort(caloHits3DFromW.begin(), caloHits3DFromW.end(),
        [&showerDirection, &showerVertex](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool
        {
            const CartesianVector positionA(pCaloHitA->GetPositionVector() - showerVertex);
            const CartesianVector positionB(pCaloHitB->GetPositionVector() - showerVertex);
            
            const float lA(showerDirection.GetDotProduct(positionA));
            const float lB(showerDirection.GetDotProduct(positionB));

            return lA < lB;
        });

    for (const CaloHit *const pHit3D : caloHits3DFromW)
    {
        const CaloHit *const pHit2D{static_cast<const CaloHit *>(pHit3D->GetParentAddress())};
        const float hitEnergy(std::fabs(pHit2D->GetElectromagneticEnergy()));
        runningEnergySum += hitEnergy;

        if ((totalEnergy > std::numeric_limits<float>::epsilon()) && 
            ((runningEnergySum / totalEnergy) > m_trueLengthEnergyFrac))
        {
            const CartesianVector displacement(pHit3D->GetPositionVector() - showerVertex);
            coreShowerLength = showerDirection.GetDotProduct(displacement);
            break;
        }
    }

    showerTreeVars.m_coreRecoLength.push_back(coreShowerLength);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetTrueLength(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticle *const pMCParticle,
     ShowerTreeVars &showerTreeVars)
{
    // Sorry for looping over matches again :( 
    CaloHitList mcHits;
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() == pMCParticle)
            mcHits = mcMatches.GetMC()->GetCaloHits();
    }

    // Get 3D positions/directions
    const CartesianVector trueShrVtx(pMCParticle->GetVertex());
    const CartesianVector trueShrDir(pMCParticle->GetMomentum().GetUnitVector());
    const CartesianVector trueShrDirSeed(trueShrVtx + (trueShrDir * 10.0f));

    // Find X% containment in longitudinal
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        float totalEnergy(0.f);
        CaloHitVector viewMCHits;
        this->GetHitsOfType(mcHits, hitType, viewMCHits, totalEnergy);

        const CartesianVector trueShrVtx2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrVtx, hitType));
        const CartesianVector trueShrDirSeed2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrDirSeed, hitType));
        const CartesianVector trueShrDir2D((trueShrDirSeed2D - trueShrVtx2D).GetUnitVector());

        // Order hits wrt l
        std::sort(viewMCHits.begin(), viewMCHits.end(),
                  [&trueShrDir2D, &trueShrVtx2D](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool
        {
            const CartesianVector positionA(pCaloHitA->GetPositionVector() - trueShrVtx2D);
            const CartesianVector positionB(pCaloHitB->GetPositionVector() - trueShrVtx2D);
            
            const float lA(trueShrDir2D.GetDotProduct(positionA));
            const float lB(trueShrDir2D.GetDotProduct(positionB));

            return lA < lB;
        });

        // Calculate 'true length'
        FloatVector &trueLength(hitType == TPC_VIEW_U ? showerTreeVars.m_coreTrueLengthFromU : 
                                hitType == TPC_VIEW_V ? showerTreeVars.m_coreTrueLengthFromV : showerTreeVars.m_coreTrueLengthFromW);

        float runningEnergySum(0.f), endpointL(-9999.f);
        for (const CaloHit *const pHit2D : viewMCHits)
        {
            const float hitEnergy(std::fabs(pHit2D->GetElectromagneticEnergy()));
            runningEnergySum += hitEnergy;

            if ((totalEnergy > std::numeric_limits<float>::epsilon()) && 
                ((runningEnergySum / totalEnergy) > m_trueLengthEnergyFrac))
            {
                const CartesianVector displacement(pHit2D->GetPositionVector() - trueShrVtx2D);
                endpointL = trueShrDir2D.GetDotProduct(displacement);
                break;
            }
        }

        if (endpointL > -9990.f)
        {
            const CartesianVector showerEndpoint2D(trueShrVtx2D + (trueShrDir2D * endpointL));
            const float scale3D((showerEndpoint2D.GetX() - trueShrVtx2D.GetX()) / trueShrDir.GetX());
            const float length3D((trueShrDir * scale3D).GetMagnitude());
            trueLength.push_back(length3D);
        }
        else
        {
            trueLength.push_back(-1.f);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetHitsOfType(const CaloHitList &inputList, const HitType hitType, CaloHitVector &outputVector, float &totalEnergy)
{
    totalEnergy = 0.f;

    for (const CaloHit *const pCaloHit : inputList)
    {
        if (pCaloHit->GetHitType() == hitType)
        {
            totalEnergy += pCaloHit->GetElectromagneticEnergy();
            outputVector.push_back(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::GetInitialRegionVars(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticle *const pMCParticle, 
    const Pfo *const pPfo, ShowerTreeVars &showerTreeVars)
{
    // Sorry for looping over matches again :( 
    CaloHitList mcHits;
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() == pMCParticle)
            mcHits = mcMatches.GetMC()->GetCaloHits();
    }

    // Get 3D positions/directions
    const CartesianVector trueShrVtx(pMCParticle->GetVertex());
    const CartesianVector trueShrDir(pMCParticle->GetMomentum().GetUnitVector());
    const CartesianVector trueShrDirSeed(trueShrVtx + (trueShrDir * m_initialRegion3D));

    // Find X% containment in longitudinal
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        float totalEnergy(0.f);
        CaloHitVector viewMCHits;
        this->GetHitsOfType(mcHits, hitType, viewMCHits, totalEnergy);
        CaloHitList viewPfoHitList;
        LArPfoHelper::GetCaloHits(pPfo, hitType, viewPfoHitList);
        CaloHitVector viewPfoHits(viewPfoHitList.begin(), viewPfoHitList.end());

        const CartesianVector trueShrVtx2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrVtx, hitType));
        const CartesianVector trueShrDirSeed2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueShrDirSeed, hitType));
        const CartesianVector trueShrDir2D((trueShrDirSeed2D - trueShrVtx2D).GetUnitVector());
        const float initialRegionL(trueShrDir2D.GetDotProduct(trueShrDirSeed2D - trueShrVtx2D));

        // Find initial region true hits
        CaloHitVector initialMCHits;
        for (const CaloHit *const pCaloHit : viewMCHits)
        {
            const float l(trueShrDir2D.GetDotProduct(pCaloHit->GetPositionVector() - trueShrVtx2D));

            if (l < initialRegionL)
                initialMCHits.push_back(pCaloHit);
        }

        // Find initial region pfo hits
        CaloHitVector initialPfoHits;
        for (const CaloHit *const pCaloHit : viewPfoHits)
        {
            const float l(trueShrDir2D.GetDotProduct(pCaloHit->GetPositionVector() - trueShrVtx2D));

            if (l < initialRegionL)
                initialPfoHits.push_back(pCaloHit);
        }

        // Get initial region completeness
        int nSharedHits(0);
        for (const CaloHit *const pCaloHit : initialMCHits)
        {
            if (std::find(initialPfoHits.begin(), initialPfoHits.end(), pCaloHit) != initialPfoHits.end())
                ++nSharedHits; 
        }

        IntVector &nInitialMCHits(hitType == TPC_VIEW_U ? showerTreeVars.m_nInitialMCHitsU : 
            hitType == TPC_VIEW_V ? showerTreeVars.m_nInitialMCHitsV : showerTreeVars.m_nInitialMCHitsW);
        IntVector &nInitialPfoHits(hitType == TPC_VIEW_U ? showerTreeVars.m_nInitialPfoHitsU : 
            hitType == TPC_VIEW_V ? showerTreeVars.m_nInitialPfoHitsV : showerTreeVars.m_nInitialPfoHitsW);
        FloatVector &completeness(hitType == TPC_VIEW_U ? showerTreeVars.m_initialCompletenessU : 
            hitType == TPC_VIEW_V ? showerTreeVars.m_initialCompletenessV : showerTreeVars.m_initialCompletenessW);
        FloatVector &purity(hitType == TPC_VIEW_U ? showerTreeVars.m_initialPurityU : 
            hitType == TPC_VIEW_V ? showerTreeVars.m_initialPurityV : showerTreeVars.m_initialPurityW);

        const float thisCompleteness(initialMCHits.size() == 0 ? -1 : static_cast<float>(nSharedHits) / initialMCHits.size());
        const float thisPurity(initialPfoHits.size() == 0 ? -1 : static_cast<float>(nSharedHits) / initialPfoHits.size());
        nInitialMCHits.push_back(initialMCHits.size());
        nInitialPfoHits.push_back(initialPfoHits.size());
        completeness.push_back(thisCompleteness);
        purity.push_back(thisPurity);
    }
}



//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::FillTree(ShowerTreeVars &showerTreeVars)
{
    FloatVector &trueShrDirX = showerTreeVars.m_trueShrDirX;
    FloatVector &trueShrDirY = showerTreeVars.m_trueShrDirY;
    FloatVector &trueShrDirZ = showerTreeVars.m_trueShrDirZ;
    FloatVector &coreTrueLengthFromU = showerTreeVars.m_coreTrueLengthFromU;
    FloatVector &coreTrueLengthFromV = showerTreeVars.m_coreTrueLengthFromV;
    FloatVector &coreTrueLengthFromW = showerTreeVars.m_coreTrueLengthFromW;
    FloatVector &recoShrVtxX = showerTreeVars.m_recoShrVtxX;
    FloatVector &recoShrVtxY = showerTreeVars.m_recoShrVtxY;
    FloatVector &recoShrVtxZ = showerTreeVars.m_recoShrVtxZ;
    FloatVector &recoShrDirX = showerTreeVars.m_recoShrDirX;
    FloatVector &recoShrDirY = showerTreeVars.m_recoShrDirY;
    FloatVector &recoShrDirZ = showerTreeVars.m_recoShrDirZ;
    FloatVector &coreRecoLength = showerTreeVars.m_coreRecoLength;
    FloatVector &recoShrLength = showerTreeVars.m_recoShrLength;
    FloatVector &recoShrDirAcc = showerTreeVars.m_recoShrDirAcc;
    FloatVector &moliereRadius = showerTreeVars.m_moliereRadius;
    IntVector &nInitialMCHitsU = showerTreeVars.m_nInitialMCHitsU;
    IntVector &nInitialMCHitsV = showerTreeVars.m_nInitialMCHitsV;
    IntVector &nInitialMCHitsW = showerTreeVars.m_nInitialMCHitsW;
    IntVector &nInitialPfoHitsU = showerTreeVars.m_nInitialPfoHitsU;
    IntVector &nInitialPfoHitsV = showerTreeVars.m_nInitialPfoHitsV;
    IntVector &nInitialPfoHitsW = showerTreeVars.m_nInitialPfoHitsW;
    FloatVector &initialCompletenessU = showerTreeVars.m_initialCompletenessU;
    FloatVector &initialCompletenessV = showerTreeVars.m_initialCompletenessV;
    FloatVector &initialCompletenessW = showerTreeVars.m_initialCompletenessW;
    FloatVector &initialPurityU = showerTreeVars.m_initialPurityU;
    FloatVector &initialPurityV = showerTreeVars.m_initialPurityV;
    FloatVector &initialPurityW = showerTreeVars.m_initialPurityW;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Run", showerTreeVars.m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Subrun", showerTreeVars.m_subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Event", showerTreeVars.m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueDirX", &trueShrDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueDirY", &trueShrDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueDirZ", &trueShrDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueCoreLengthFromU", &coreTrueLengthFromU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueCoreLengthFromV", &coreTrueLengthFromV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueCoreLengthFromW", &coreTrueLengthFromW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoVtxX", &recoShrVtxX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoVtxY", &recoShrVtxY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoVtxZ", &recoShrVtxZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoDirX", &recoShrDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoDirY", &recoShrDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoDirZ", &recoShrDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoCoreLength", &coreRecoLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoLength", &recoShrLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_DirAcc", &recoShrDirAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_MoliereRadius", &moliereRadius));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialMCHitsU", &nInitialMCHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialMCHitsV", &nInitialMCHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialMCHitsW", &nInitialMCHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialPfoHitsU", &nInitialPfoHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialPfoHitsV", &nInitialPfoHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialPfoHitsW", &nInitialPfoHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialCompletenessU", &initialCompletenessU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialCompletenessV", &initialCompletenessV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialCompletenessW", &initialCompletenessW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialPurityU", &initialPurityU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialPurityV", &initialPurityV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_InitialPurityW", &initialPurityW));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ShowerTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InitialRegion", m_initialRegion3D));


    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrueLengthEnergyFrac", m_trueLengthEnergyFrac));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
