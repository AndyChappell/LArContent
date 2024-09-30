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
#include "larpandoracontent/LArObjects/LArGraph.h"

#include "larpandoradlcontent/LArVertex/DlSageVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlSageVertexingAlgorithm::DlSageVertexingAlgorithm() :
    m_trainingMode{false},
    m_event{-1},
    m_nClasses{0},
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
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        LArGraph graph;
        graph.MakeGraph(*pCaloHitList);
        const LArGraph::EdgeVector &edgeVector{graph.GetEdges()};
        int id{0};
        std::map<const CaloHit *, int> nodeIdMap;
        IntVector node0Vector, node1Vector;
        float x{0.f}, z{0.f}, adc{0.f};
        for (const LArGraph::Edge *pEdge : edgeVector)
        {
            const CaloHit *pNode0{pEdge->m_v0};
            const CaloHit *pNode1{pEdge->m_v1};
            if (nodeIdMap.find(pNode0) == nodeIdMap.end())
            {
                nodeIdMap[pNode0] = id++;
                // Update and fill the node tree - think about scaling here, probably want to get to a final feature set
                // at this stage, rather than computing things in Python
                x = pNode0->GetPositionVector().GetX();
                z = pNode0->GetPositionVector().GetZ();
                adc = pNode0->GetInputEnergy();
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "idx", id));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "x", x));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "z", z));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "adc", adc));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_nodeTreeName));
            }
            if (nodeIdMap.find(pNode1) == nodeIdMap.end())
            {
                nodeIdMap[pNode1] = id++;
                // Update and fill the node tree
                x = pNode1->GetPositionVector().GetX();
                z = pNode1->GetPositionVector().GetZ();
                adc = pNode1->GetInputEnergy();
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "idx", id));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "x", x));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "z", z));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_nodeTreeName, "adc", adc));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_nodeTreeName));
            }
            // Update the edge tree
            node0Vector.emplace_back(nodeIdMap.at(pNode0));
            node1Vector.emplace_back(nodeIdMap.at(pNode1));
        }
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

