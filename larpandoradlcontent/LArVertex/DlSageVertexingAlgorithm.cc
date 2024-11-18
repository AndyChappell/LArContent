/**
 *  @file   larpandoradlcontent/LArVertex/DlSageVertexingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning vertexing algorithm.
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoradlcontent/LArVertex/DlSageVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlSageVertexingAlgorithm::DlSageVertexingAlgorithm() :
    m_trainingMode{false},
    m_event{-1},
    m_nClasses{0},
    m_nEdges{5},
    m_maxSecondaryDistance{3.f},
    m_maxSecondaryCosine{0.996f},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

DlSageVertexingAlgorithm::~DlSageVertexingAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_edgeTreeName, m_graphFileName, "RECREATE"));
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_nodeTreeName, m_graphFileName, "UPDATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "DlSageVertexingAlgorithm: Failed to save ROOT tree" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlSageVertexingAlgorithm::Visualize(const LArGraph &graph) const
{
#ifdef MONITORING
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    for (const auto pEdge : graph.GetEdges())
    {
        const CartesianVector &p0{pEdge->m_v0->GetPositionVector()}, &p1{pEdge->m_v1->GetPositionVector()};
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &p0, &p1, "edge", BLUE, 1, 1));
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSageVertexingAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlSageVertexingAlgorithm::PrepareTrainingSample()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    if (!pMCParticleList || pMCParticleList->empty())
        return STATUS_CODE_NOT_FOUND;

    const MCParticle *pNeutrino{nullptr};
    for (const MCParticle *const pMC : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMC) && pMC->GetParentList().empty())
            pNeutrino = pMC;
    }
    if (!pNeutrino)
        return STATUS_CODE_SUCCESS;

    const CartesianVector &vertex(pNeutrino->GetVertex());
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    const float xVtx{vertex.GetX()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        float zVtx{0.f};
        switch (pCaloHitList->front()->GetHitType())
        {
            case TPC_VIEW_U:
                zVtx = transform->YZtoU(vertex.GetY(), vertex.GetZ());
                break;
            case TPC_VIEW_V:
                zVtx = transform->YZtoV(vertex.GetY(), vertex.GetZ());
                break;
            case TPC_VIEW_W:
                zVtx = transform->YZtoW(vertex.GetY(), vertex.GetZ());
                break;
            default:
                continue;
        }

        LArGraph graph(false, this->m_nEdges, this->m_maxSecondaryCosine, this->m_maxSecondaryDistance);
        graph.MakeGraph(*pCaloHitList);
        //this->Visualize(graph);
        const LArGraph::EdgeVector &edgeVector{graph.GetEdges()};
        int id{0};
        std::map<const CaloHit *, int> nodeIdMap;
        IntVector node0Vector, node1Vector, nodeIdVector;
        FloatVector xVector, zVector, adcVector, distanceVector;
        for (const LArGraph::Edge *pEdge : edgeVector)
        {
            const CaloHit *pNode0{pEdge->m_v0};
            const CaloHit *pNode1{pEdge->m_v1};
            for (const CaloHit *pNode : {pNode0, pNode1})
            {
                if (nodeIdMap.find(pNode) == nodeIdMap.end())
                {
                    nodeIdMap[pNode] = id;
                    nodeIdVector.emplace_back(id);
                    xVector.emplace_back(pNode->GetPositionVector().GetX());
                    zVector.emplace_back(pNode->GetPositionVector().GetZ());
                    adcVector.emplace_back(pNode->GetInputEnergy());
                    const float dx{pNode->GetPositionVector().GetX() - xVtx};
                    const float dz{pNode->GetPositionVector().GetZ() - zVtx};
                    const float dr{std::sqrt(dx * dx + dz * dz)};
                    distanceVector.emplace_back(dr);
                    ++id;
                }
            }
            // Update the edge tree
            node0Vector.emplace_back(nodeIdMap.at(pNode0));
            node1Vector.emplace_back(nodeIdMap.at(pNode1));
            node0Vector.emplace_back(nodeIdMap.at(pNode1));
            node1Vector.emplace_back(nodeIdMap.at(pNode0));
        }
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "id_vector", &nodeIdVector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "x_vector", &xVector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "z_vector", &zVector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "adc_vector", &adcVector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "distance_vector", &distanceVector));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_nodeTreeName));

        // Provisionally this will need to be done for each of u, v, and w, so we'll need unique tree names
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_edgeTreeName, "node0_vector", &node0Vector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_edgeTreeName, "node1_vector", &node1Vector));
        // Might want to add edge attributes here - like distance (factoring in direction)
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_edgeTreeName));
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSageVertexingAlgorithm::Infer()
{
    ++m_event;

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSageVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DistanceThresholds", m_thresholds));
    m_nClasses = m_thresholds.size() - 1;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxNodeEdges", m_nEdges));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSecondaryDistance", m_maxSecondaryDistance));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSecondaryCosine", m_maxSecondaryCosine));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "GraphFileName", m_graphFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "EdgeTreeName", m_edgeTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NodeTreeName", m_nodeTreeName));
    }
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameU", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelU);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameV", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelV);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameW", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelW);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content

