/**
 *  @file   larpandoradlcontent/LArVertex/DlSparseVertexingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning vertexing algorithm.
 *
 *  $Log: $
 */

#include <chrono>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoradlcontent/LArVertex/DlSparseVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlSparseVertexingAlgorithm::DlSparseVertexingAlgorithm() :
    m_caloHitListName{"CaloHitList2D"},
    m_rootFileName{"sparse_vertexing.root"},
    m_rootTreeName{"sparse_vertexing"},
    m_trainingMode{false}
{
}

DlSparseVertexingAlgorithm::~DlSparseVertexingAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlSparseVertexingAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSparseVertexingAlgorithm::Run()
{
    if (!m_trainingMode)
        return this->Infer();
    else
        return this->PrepareTrainingSample();

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSparseVertexingAlgorithm::PrepareTrainingSample()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    // Quiet exit if the list is absent or empty
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    CartesianVector trueVertex3D(0, 0, 0);
    // Quiet exit if no true vertex could be found
    if (!LArMCParticleHelper::GetTrueVertex(pMCParticleList, trueVertex3D))
        return STATUS_CODE_SUCCESS;

    const HitType view{pCaloHitList->front()->GetHitType()};
    const float pitch{LArGeometryHelper::GetWirePitch(this->GetPandora(), view)};
    if (pitch <= 0.f)
    {
        std::cout << "Error - DlSparseVertexingAlgorithm: Unable to get wire pitch for view " << view << std::endl;
        return STATUS_CODE_NOT_ALLOWED;
    }
    // Zero initialisation is fine, the function will overwrite it
    float xMin{0}, xMax{0}, zMin{0}, zMax{0};
    this->GetHitRegion(*pCaloHitList, xMin, xMax, zMin, zMax, std::numeric_limits<float>::epsilon());
    //const float xIndexMin{0.f}, xIndexMax{std::floor((xMax - xMin) / pitch)};
    //const float zIndexMin{0.f}, zIndexMax{std::floor((zMax - zMin) / pitch)};
    //const float scale{std::sqrt((xIndexMax - xIndexMin) * (xIndexMax - xIndexMin) + (zIndexMax - zIndexMin) * (zIndexMax - zIndexMin))};

    // Get the true vertex location and convert to pixel coordinates
    float vx{0.f}, vu{0.f}, vv{0.f}, vw{0.f};
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    LArVertexHelper::GetTrueVertexPosition(trueVertex3D, transform, vx, vu, vv, vw);
    float vz{0.f};
    switch (view)
    {
        case TPC_VIEW_U:
            vz = vu;
            break;
        case TPC_VIEW_V:
            vz = vv;
            break;
        case TPC_VIEW_W:
            vz = vw;
            break;
        default:
            std::cout << "Error - DlSparseVertexingAlgorithm: Unknown view " << view << std::endl;
            return STATUS_CODE_NOT_ALLOWED;
    }

    const float vxIndex{std::floor((vx - xMin) / pitch)};
    const float vzIndex{std::floor((vz - zMin) / pitch)};

    std::map<std::pair<int, int>, std::tuple<float, float, int>> featureMap;
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
        const float width{pCaloHit->GetCellSize1() / pitch};
        const int xIndex{static_cast<int>(std::floor((x - xMin) / pitch))};
        const int zIndex{static_cast<int>(std::floor((z - zMin) / pitch))};
        if (featureMap.find({xIndex, zIndex}) == featureMap.end())
        {
            //float dr{std::sqrt((vxIndex - xIndex) * (vxIndex - xIndex) + (vzIndex - zIndex) * (vzIndex - zIndex)) / scale};
            float dr{std::sqrt((vxIndex - xIndex) * (vxIndex - xIndex) + (vzIndex - zIndex) * (vzIndex - zIndex))};
            int cls{this->GetClassFromDistance(dr)};
            featureMap[{xIndex, zIndex}] = std::make_tuple(pCaloHit->GetMipEquivalentEnergy(), width, cls);
        }
        else
        {
            // Pick the maximum energy deposit in the pixel
            if (std::get<0>(featureMap[{xIndex, zIndex}]) < pCaloHit->GetMipEquivalentEnergy())
                featureMap[{xIndex, zIndex}] = std::make_tuple(pCaloHit->GetMipEquivalentEnergy(), width, std::get<2>(featureMap[{xIndex, zIndex}]));
        }
    }
    // Consolidate coordinates and features
    IntVector xx, zz;
    FloatVector adc, width, distance;
    for (const auto &[key, value] : featureMap)
    {
        xx.emplace_back(key.first);
        zz.emplace_back(key.second);
        adc.emplace_back(std::get<0>(value));
        width.emplace_back(std::get<1>(value));
        distance.emplace_back(std::get<2>(value));
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "xx", &xx));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "zz", &zz));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "adc", &adc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "width", &width));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "distance", &distance));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSparseVertexingAlgorithm::Infer()
{
    LArDLHelper::TorchModel model;
    LArDLHelper::LoadModel(LArFileHelper::FindFileInPath("PandoraNetworkData/sparse_unet.pt", "FW_SEARCH_PATH"), model);
    model.eval();

    torch::Tensor coords = torch::randint(0, 255, {600, 2}, torch::kInt64);
    torch::Tensor features = torch::randn({600, 3}, torch::kFloat32);

    auto start = std::chrono::high_resolution_clock::now();

    LArDLHelper::TorchInputVector inputs;
    inputs.push_back(coords);
    inputs.push_back(features);

    auto output = model.forward(inputs);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::micro> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " microseconds" << std::endl;

    std::cout << "Inferring vertex location" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSparseVertexingAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax, const float padding) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        xMin = std::min(pos.GetX(), xMin);
        xMax = std::max(pos.GetX(), xMax);
        zMin = std::min(pos.GetZ(), zMin);
        zMax = std::max(pos.GetZ(), zMax);
    }

    if (padding > 0.f)
    {
        xMin -= padding;
        xMax += padding;
        zMin -= padding;
        zMax += padding;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

int DlSparseVertexingAlgorithm::GetClassFromDistance(float distance) const
{
    int left{0}, right{static_cast<int>(m_thresholds.size()) - 1};

    while (left <= right)
    {
        if (left == right)
            return left;
        const int mid{left + (right - left) / 2};
        if (m_thresholds[mid] == distance)
            return mid;
        else if (m_thresholds[mid] < distance)
        {
            if (left < mid)
                left = mid;
            else
                right = mid;
        }
        else
            right = mid;
    }

    return left;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSparseVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (!m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "Thresholds", m_thresholds));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
