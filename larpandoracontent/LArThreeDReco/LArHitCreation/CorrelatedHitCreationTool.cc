/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/CorrelatedHitCreationTool.cc
 *
 *  @brief  Implementation of the hit creation base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/CorrelatedHitCreationTool.h"

using namespace pandora;

namespace lar_content
{

CorrelatedHitCreationTool::CorrelatedHitCreationTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CorrelatedHitCreationTool::~CorrelatedHitCreationTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CorrelatedHitCreationTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo,
    const CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
    (void)pPfo;

    std::map<const CaloHit *, bool> hitMap;
    for (const CaloHit *const pCaloHit : inputTwoDHits)
        hitMap[pCaloHit] = true;

    CaloHitList referenceCaloHits;
    for (const CaloHit *const pCaloHit : inputTwoDHits)
    {
        const LArCaloHit *const pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        if (!pLArCaloHit)
            continue;
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            // A single sibling in the same PFO allows us to create a 3D hit for the current hit
            const CaloHit *const pSiblingHit{pLArCaloHit->GetSiblingHit(view)};
            if (pSiblingHit && hitMap.find(pSiblingHit) != hitMap.end())
            {
                referenceCaloHits.emplace_back(pLArCaloHit);
                break;
            }
        }
    }
    std::cout << "Ref hits: " << referenceCaloHits.size() << std::endl;

    for (const CaloHit *const pRefCaloHit : referenceCaloHits)
    {
        const LArCaloHit *const pLArRefCaloHit{dynamic_cast<const LArCaloHit *>(pRefCaloHit)};
        if (!pLArRefCaloHit)
            continue;
        ProtoHit protoHit{pLArRefCaloHit};
        const LArCaloHit *pSiblingHit1{nullptr}, *pSiblingHit2{nullptr};
        switch (pLArRefCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                pSiblingHit1 = dynamic_cast<const LArCaloHit *>(pLArRefCaloHit->GetSiblingHit(TPC_VIEW_V));
                pSiblingHit2 = dynamic_cast<const LArCaloHit *>(pLArRefCaloHit->GetSiblingHit(TPC_VIEW_W));
                break;
            case TPC_VIEW_V:
                pSiblingHit1 = dynamic_cast<const LArCaloHit *>(pLArRefCaloHit->GetSiblingHit(TPC_VIEW_U));
                pSiblingHit2 = dynamic_cast<const LArCaloHit *>(pLArRefCaloHit->GetSiblingHit(TPC_VIEW_W));
                break;
            case TPC_VIEW_W:
                pSiblingHit1 = dynamic_cast<const LArCaloHit *>(pLArRefCaloHit->GetSiblingHit(TPC_VIEW_U));
                pSiblingHit2 = dynamic_cast<const LArCaloHit *>(pLArRefCaloHit->GetSiblingHit(TPC_VIEW_V));
                break;
            default:
                break;
        }

        if (pSiblingHit1 && pSiblingHit2)
        {
            const CartesianVector &pos0{pLArRefCaloHit->GetPositionVector()};
            const CartesianVector &pos1{pSiblingHit1->GetPositionVector()};
            const CartesianVector &pos2{pSiblingHit2->GetPositionVector()};
            CartesianVector pos3D(0, 0, 0);
            float chi2{0.f};
            // Currently slightly different to the Correlated approach due to lack of error consideration
            LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), pLArRefCaloHit->GetHitType(), pSiblingHit1->GetHitType(),
                pSiblingHit2->GetHitType(), pos0, pos1, pos2, pos3D, chi2);
            protoHit.SetPosition3D(pos3D, chi2);
            protoHitVector.emplace_back(protoHit);
        }
        else if (pSiblingHit1)
        {
            const CartesianVector &pos0{pLArRefCaloHit->GetPositionVector()};
            const CartesianVector &pos1{pSiblingHit1->GetPositionVector()};
            CartesianVector pos3D(0, 0, 0);
            float chi2{0.f};
            // Currently slightly different to the Correlated approach due to lack of error consideration
            LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), pLArRefCaloHit->GetHitType(), pSiblingHit1->GetHitType(),
                pos0, pos1, pos3D, chi2);
            protoHit.SetPosition3D(pos3D, chi2);
            protoHitVector.emplace_back(protoHit);
        }
        else if (pSiblingHit2)
        {
            const CartesianVector &pos0{pLArRefCaloHit->GetPositionVector()};
            const CartesianVector &pos2{pSiblingHit2->GetPositionVector()};
            CartesianVector pos3D(0, 0, 0);
            float chi2{0.f};
            // Currently slightly different to the Correlated approach due to lack of error consideration
            LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), pLArRefCaloHit->GetHitType(), pSiblingHit2->GetHitType(),
                pos0, pos2, pos3D, chi2);
            protoHit.SetPosition3D(pos3D, chi2);
            protoHitVector.emplace_back(protoHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CorrelatedHitCreationTool::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    return HitCreationBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
