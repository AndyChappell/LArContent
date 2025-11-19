/**
 *  @file   larpandoracontent/LArMetrics/ShowerValidationTool.cc
 *
 *  @brief  Implementation of the shower validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMetrics/ShowerValidationTool.h"

using namespace pandora;

namespace lar_content
{

ShowerValidationTool::ShowerValidationTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticleVector &targetMC, const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ShowerTreeVars showerTreeVars;

    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        const MCParticle *const pMC(targetMC.at(i));
        const CartesianVector trueDir(pMC->GetMomentum().GetUnitVector());

        showerTreeVars.m_trueShrDirX.push_back(trueDir.GetX());
        showerTreeVars.m_trueShrDirY.push_back(trueDir.GetY());
        showerTreeVars.m_trueShrDirZ.push_back(trueDir.GetZ());

        const Pfo *const pPfo(bestRecoMatch.at(i));

        if (pPfo)
        {
            CartesianVector recoShrVtx(0.f, 0.f, 0.f), recoShrDir(0.f, 0.f, 0.f);
            this->GetShowerDirection(pPfo, recoShrVtx, recoShrDir);

            showerTreeVars.m_recoShrVtxX.push_back(recoShrVtx.GetX());
            showerTreeVars.m_recoShrVtxY.push_back(recoShrVtx.GetY());
            showerTreeVars.m_recoShrVtxZ.push_back(recoShrVtx.GetZ());
            showerTreeVars.m_recoShrDirX.push_back(recoShrDir.GetX());
            showerTreeVars.m_recoShrDirY.push_back(recoShrDir.GetY());
            showerTreeVars.m_recoShrDirZ.push_back(recoShrDir.GetZ());
            showerTreeVars.m_recoShrDirAcc.push_back(trueDir.GetOpeningAngle(recoShrDir));

            this->GetMoliere(pPfo, recoShrVtx, recoShrDir, showerTreeVars);
        }
        else
        {
            showerTreeVars.m_recoShrVtxX.push_back(-9999.f);
            showerTreeVars.m_recoShrVtxY.push_back(-9999.f);
            showerTreeVars.m_recoShrVtxZ.push_back(-9999.f);
            showerTreeVars.m_recoShrDirX.push_back(-9999.f);
            showerTreeVars.m_recoShrDirY.push_back(-9999.f);
            showerTreeVars.m_recoShrDirZ.push_back(-9999.f);
            showerTreeVars.m_recoShrDirAcc.push_back(-1.f);
            showerTreeVars.m_moliereRadius.push_back(-1.f);
        }
    }

    this->FillTree(showerTreeVars);
}

//------------------------------------------------------------------------------------------------------------------------------------------

// Replicate PandoraShower fitting
void ShowerValidationTool::GetShowerDirection(const Pfo *const pPfo, CartesianVector &showerVertex, CartesianVector &showerDirection)
{
    CartesianPointVector positions3D;
    LArPfoHelper::GetCoordinateVector(pPfo, TPC_3D, positions3D);
    const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(pPfo));
    const CartesianVector vertexPosition(pRecoVertex->GetPosition());
    // ATTN: do they always have a vertex?
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

    std::sort(caloHits3DFromW.begin(), caloHits3DFromW.end(),
        [&showerDirection, &showerVertex](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool
        {
            const CartesianVector positionA(pCaloHitA->GetPositionVector() - showerVertex);
            const CartesianVector positionB(pCaloHitB->GetPositionVector() - showerVertex);
            
            const float tA(showerDirection.GetCrossProduct(positionA).GetMagnitude());
            const float tB(showerDirection.GetCrossProduct(positionB).GetMagnitude());

            return tA < tB;
        });

    // Now do Molliere
    float runningEnergySum(0.f), moliereRadius(-1.f);

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerValidationTool::FillTree(ShowerTreeVars &showerTreeVars)
{
    FloatVector &trueShrDirX = showerTreeVars.m_trueShrDirX;
    FloatVector &trueShrDirY = showerTreeVars.m_trueShrDirY;
    FloatVector &trueShrDirZ = showerTreeVars.m_trueShrDirZ;
    FloatVector &recoShrVtxX = showerTreeVars.m_recoShrVtxX;
    FloatVector &recoShrVtxY = showerTreeVars.m_recoShrVtxY;
    FloatVector &recoShrVtxZ = showerTreeVars.m_recoShrVtxZ;
    FloatVector &recoShrDirX = showerTreeVars.m_recoShrDirX;
    FloatVector &recoShrDirY = showerTreeVars.m_recoShrDirY;
    FloatVector &recoShrDirZ = showerTreeVars.m_recoShrDirZ;
    FloatVector &recoShrDirAcc = showerTreeVars.m_recoShrDirAcc;
    FloatVector &moliereRadius = showerTreeVars.m_moliereRadius;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueDirX", &trueShrDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueDirY", &trueShrDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_TrueDirZ", &trueShrDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoVtxX", &recoShrVtxX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoVtxY", &recoShrVtxY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoVtxZ", &recoShrVtxZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoDirX", &recoShrDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoDirY", &recoShrDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_RecoDirZ", &recoShrDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_DirAcc", &recoShrDirAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ShowerTree", "Shr_MoliereRadius", &moliereRadius));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ShowerTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerValidationTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
