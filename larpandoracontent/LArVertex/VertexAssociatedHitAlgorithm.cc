/**
 *  @file   larpandoracontent/LArTwoDReco/LArAssociatedHit/VertexAssociatedHitAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArVertex/VertexAssociatedHitAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

VertexAssociatedHitAlgorithm::VertexAssociatedHitAlgorithm() :
    m_caloHitListName{""},
    m_vertexListName{""}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssociatedHitAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const VertexList *pVertexList{nullptr};
    StatusCode status{PandoraContentApi::GetList(*this, m_vertexListName, pVertexList)};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, status);

    if (status == STATUS_CODE_NOT_INITIALIZED)
        return STATUS_CODE_SUCCESS;

    if (!pCaloHitList || pCaloHitList->empty() || !pVertexList || pVertexList->empty())
        return STATUS_CODE_SUCCESS;

    this->IdentifyAssociatedHits(*pCaloHitList, *pVertexList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexAssociatedHitAlgorithm::IdentifyAssociatedHits(const CaloHitList &caloHitList, const VertexList &vertexList) const
{
    HitType view{caloHitList.front()->GetHitType()};
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};

    // Gather the hits within a given radii of the vertices
    std::set<const CaloHit *> vetoedHitSet;
    for (const Vertex *const pVertex : vertexList)
    {
        const CartesianVector &pos3D{pVertex->GetPosition()};
        CartesianVector pos(0, 0, 0);
        switch (view)
        {
            case TPC_VIEW_U:
                pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoU(pos3D.GetY(), pos3D.GetZ())));
                break;
            case TPC_VIEW_V:
                pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoV(pos3D.GetY(), pos3D.GetZ())));
                break;
            case TPC_VIEW_W:
                pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoW(pos3D.GetY(), pos3D.GetZ())));
                break;
            default:
                return;
        }

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
            const float distanceSquared{hitPos.GetDistanceSquared(pos)};
            if (distanceSquared <= 25.f)
                vetoedHitSet.insert(pCaloHit);
        }
    }

    CaloHitList vetoedHits;
    for (const CaloHit *const pCaloHit :vetoedHitSet)
    {
        vetoedHits.emplace_back(pCaloHit);
    }

    CaloHitList retainedHits;
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (vetoedHitSet.find(pCaloHit) == vetoedHitSet.end())
        {
            retainedHits.emplace_back(pCaloHit);
        }
    }

    PandoraContentApi::SaveList(*this, vetoedHits, m_vetoedHitListName);
    PandoraContentApi::SaveList(*this, retainedHits, m_retainedHitListName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssociatedHitAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VetoedHitListName", m_vetoedHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RetainedHitListName", m_retainedHitListName));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
