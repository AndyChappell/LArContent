/**
 *  @file   larpandoracontent/LArCheating/CheatingSecondaryVertexCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArCheating/CheatingSecondaryVertexCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingSecondaryVertexCreationAlgorithm::CheatingSecondaryVertexCreationAlgorithm() : m_replaceCurrentVertexList(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingSecondaryVertexCreationAlgorithm::Run()
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));
    CartesianPointVector vertices;
    VisibleParticleMap visibleParticleMap;
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMatchedMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            visibleParticleMap[pMatchedMC] = true;
        }
        catch (...)
        {
        }
    }

    for (const MCParticle *mc : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(mc))
        {
            const MCParticleList &primaries{mc->GetDaughterList()};
            for (const MCParticle *pPrimary : primaries)
            {
                if (visibleParticleMap.find(pPrimary) == visibleParticleMap.end())
                    continue;

                CartesianVector vertex(0, 0, 0);
                if (this->GetClearSecondaryVertex(pPrimary, visibleParticleMap, vertex))
                    vertices.emplace_back(vertex);
            }
        }
    }
    if (!vertices.empty())
    {
        // Make vertices
        const VertexList *pVertexList(nullptr);
        std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        for (const CartesianVector &vertex : vertices)
        {
            //std::cout << "Found secondary at " << vertex << std::endl;
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = vertex;
            parameters.m_vertexLabel = VERTEX_INTERACTION;
            parameters.m_vertexType = VERTEX_3D;

            const Vertex *pVertex(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
        }

        if (!pVertexList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool CheatingSecondaryVertexCreationAlgorithm::GetClearSecondaryVertex(const MCParticle *const pParent, const VisibleParticleMap &visibleParticleMap,
    CartesianVector &vertex) const
{
    const MCParticleList &secondaries{pParent->GetDaughterList()};
    const MCParticle *pLastSecondary{nullptr};
    int numVisible{0};
    for (const MCParticle *pSecondary : secondaries)
    {
        if (visibleParticleMap.find(pSecondary) != visibleParticleMap.end())
        {
            pLastSecondary = pSecondary;
            ++numVisible;
            break;
        }
    }
    if (numVisible > 0)
    {
        if (numVisible == 1 && pParent->GetParticleId() == pLastSecondary->GetParticleId())
        {
            const CartesianVector &a{pParent->GetMomentum().GetUnitVector()};
            const CartesianVector &b{pLastSecondary->GetMomentum().GetUnitVector()};
            const float costheta{a.GetDotProduct(b)};
            if (costheta <= 0.985f)
            {
                const CartesianVector &endPoint{pParent->GetEndpoint()};
                vertex.SetValues(endPoint.GetX(), endPoint.GetY(), endPoint.GetZ());

                return true;
            }
            else
            {
                return this->GetClearSecondaryVertex(pLastSecondary, visibleParticleMap, vertex);
            }
        }
        else
        {
            const CartesianVector &endPoint{pParent->GetEndpoint()};
            vertex.SetValues(endPoint.GetX(), endPoint.GetY(), endPoint.GetZ());

            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingSecondaryVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
