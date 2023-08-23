/**
 *  @file   larpandoracontent/LArMonitoring/SecondaryVertexMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the secondary vertex monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/SecondaryVertexMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoracontent/LArObjects/LArEventTopology.h"

using namespace pandora;

namespace lar_content
{

SecondaryVertexMonitoringAlgorithm::SecondaryVertexMonitoringAlgorithm() :
    m_writeFile{false},
    m_caloHitListName{"CaloHitList2D"},
    m_vertexListName{"NeutrinoVertices3D"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SecondaryVertexMonitoringAlgorithm::~SecondaryVertexMonitoringAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryVertexMonitoringAlgorithm::Run()
{
    this->AssessVertices();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryVertexMonitoringAlgorithm::AssessVertices() const
{
#ifdef MONITORING
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList2D));
    const VertexList *pVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    LArEventTopology eventTopology(*pCaloHitList2D);
    eventTopology.ConstructVisibleHierarchy();
    eventTopology.PruneHierarchy();
    CartesianPointVector trueVertices;
    eventTopology.GetVertices(trueVertices);

    std::map<const CartesianVector *, CartesianPointVector> vertexMap;
    for (const CartesianVector &trueVertex : trueVertices)
    {
        float bestDistance{std::numeric_limits<float>::max()};
        float bestDx{bestDistance}, bestDy{bestDistance}, bestDz{bestDistance};
        for (const Vertex *pRecoVertex : *pVertexList)
        {
            const CartesianVector &recoPosition{pRecoVertex->GetPosition()};
            const CartesianVector &displacement{recoPosition - trueVertex};
            const float dx{displacement.GetX()}, dy{displacement.GetY()}, dz{displacement.GetZ()};
            const float distance{displacement.GetMagnitude()};
            if (distance < bestDistance && distance < 5.f)
            {
                bestDistance = distance;
                bestDx = dx;
                bestDy = dy;
                bestDz = dz;
                vertexMap[&trueVertex].insert(vertexMap[&trueVertex].begin(), recoPosition);
            }
            else if (distance < 5.f)
            {
                vertexMap[&trueVertex].emplace_back(recoPosition);
            }
        }
        const int nMatches{static_cast<int>(vertexMap[&trueVertex].size())};
        if (nMatches == 0)
            bestDistance = -1;
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", bestDx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", bestDy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", bestDz));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", bestDistance));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
    }
#endif

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryVertexMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
