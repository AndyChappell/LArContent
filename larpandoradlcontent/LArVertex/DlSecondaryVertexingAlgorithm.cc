/**
 *  @file   larpandoradlcontent/LArVertex/DlSecondaryVertexingAlgorithm.cc
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
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoracontent/LArObjects/LArEventTopology.h"

#include "larpandoradlcontent/LArVertex/DlSecondaryVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlSecondaryVertexingAlgorithm::DlSecondaryVertexingAlgorithm() :
    m_trainingMode{false},
    m_trainingOutputFile{""},
    m_event{-1},
    m_pass{1},
    m_nClasses{0},
    m_height{256},
    m_width{256},
    m_driftStep{0.5f},
    m_visualise{false},
    m_writeTree{false},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

DlSecondaryVertexingAlgorithm::~DlSecondaryVertexingAlgorithm()
{
    if (m_writeTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "VertexAssessmentAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::Run()
{
    ++m_event;

    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlSecondaryVertexingAlgorithm::PrepareTrainingSample()
{
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    LArEventTopology eventTopology(*pCaloHitList2D);
    eventTopology.ConstructVisibleHierarchy();
    eventTopology.PruneHierarchy();
    CartesianPointVector vertices;
    eventTopology.GetVertices(vertices);

    // Only train on events where there is a vertex in the fiducial volume
    bool hasFiducialVertex{false};
    for (const CartesianVector &vertex : vertices)
    {
        if (LArVertexHelper::IsInFiducialVolume(this->GetPandora(), vertex, "dune_fd_hd"))
        {
            hasFiducialVertex = true;
            break;
        }
    }

    if (!hasFiducialVertex)
        return STATUS_CODE_SUCCESS;

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMCToHitsMap(mcToHitsMap));
    MCParticleList hierarchy;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CompleteMCHierarchy(mcToHitsMap, hierarchy));

    // Get boundaries for hits and make x dimension common
    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};

        const std::string trainingFilename{m_trainingOutputFile + "_" + listname + ".csv"};
        const unsigned long nVertices{vertices.size()};
        unsigned long nHits{0};

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(m_event));
        featureVector.emplace_back(static_cast<double>(nVertices));
        // Vertices
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        for (const CartesianVector &vertex : vertices)
        {
            switch (view)
            {
                case TPC_VIEW_U:
                    featureVector.emplace_back(vertex.GetX());
                    featureVector.emplace_back(transform->YZtoU(vertex.GetY(), vertex.GetZ()));
                    break;
                case TPC_VIEW_V:
                    featureVector.emplace_back(vertex.GetX());
                    featureVector.emplace_back(transform->YZtoV(vertex.GetY(), vertex.GetZ()));
                    break;
                default:
                    featureVector.emplace_back(vertex.GetX());
                    featureVector.emplace_back(transform->YZtoW(vertex.GetY(), vertex.GetZ()));
                    break;
            }
        }

        // Retain the hit region
        featureVector.emplace_back(driftMin);
        featureVector.emplace_back(driftMax);
        featureVector.emplace_back(wireMin[view]);
        featureVector.emplace_back(wireMax[view]);

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            featureVector.emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
            featureVector.emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
            featureVector.emplace_back(static_cast<double>(pCaloHit->GetMipEquivalentEnergy()));
            ++nHits;
        }
        featureVector.insert(featureVector.begin() + 2 + 2 * nVertices + 4, static_cast<double>(nHits));
        // Only write out the feature vector if there were enough hits in the region of interest
        if (nHits > 10)
            LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::Infer()
{
    // Get boundaries for hits and make x dimension common
    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    CanvasViewMap canvases;
    canvases[TPC_VIEW_U] = nullptr;
    canvases[TPC_VIEW_V] = nullptr;
    canvases[TPC_VIEW_W] = nullptr;
    CartesianPointVector vertexCandidatesU, vertexCandidatesV, vertexCandidatesW;
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        LArDLHelper::TorchInput input;
        PixelVector pixelVector;
        this->MakeNetworkInputFromHits(*pCaloHitList, view, driftMin, driftMax, wireMin[view], wireMax[view], input, pixelVector);

        // Run the input through the trained model
        LArDLHelper::TorchInputVector inputs;
        inputs.push_back(input);
        LArDLHelper::TorchOutput output;
        switch (view)
        {
            case TPC_VIEW_U:
                LArDLHelper::Forward(m_modelU, inputs, output);
                break;
            case TPC_VIEW_V:
                LArDLHelper::Forward(m_modelV, inputs, output);
                break;
            default:
                LArDLHelper::Forward(m_modelW, inputs, output);
                break;
        }

        int colOffset{0}, rowOffset{0}, canvasWidth{m_width}, canvasHeight{m_height};
        this->GetCanvasParameters(output, pixelVector, colOffset, rowOffset, canvasWidth, canvasHeight);
        canvases[view] = new Canvas(view, canvasWidth, canvasHeight, colOffset, rowOffset, driftMin, driftMax, wireMin[view], wireMax[view]);
        // we want the maximum value in the num_classes dimension (1) for every pixel
        auto classes{torch::argmax(output, 1)};
        // the argmax result is a 1 x height x width tensor where each element is a class id
        auto classesAccessor{classes.accessor<long, 3>()};
        const double scaleFactor{std::sqrt(m_height * m_height + m_width * m_width)};
        std::map<int, bool> haveSeenMap;
        for (const auto &[row, col] : pixelVector)
        {
            const auto cls{classesAccessor[0][row][col]};
            if (cls > 0 && cls < m_nClasses)
            {
                const int inner{static_cast<int>(std::round(std::ceil(scaleFactor * m_thresholds[cls - 1])))};
                const int outer{static_cast<int>(std::round(std::ceil(scaleFactor * m_thresholds[cls])))};
                this->DrawRing(canvases[view]->m_canvas, row + rowOffset, col + colOffset, inner, outer, 1.f / (outer * outer - inner * inner));
            }
        }
    }

    CartesianPointVector vertexVector;
    this->GetNetworkVertices(canvases, vertexVector);

    if (!vertexVector.empty())
    {
        StatusCode status{this->MakeCandidateVertexList(vertexVector)};

        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            if (canvases[view])
                delete canvases[view];
        }

        return status;
    }
    else
    {
        std::cout << "Insufficient 2D vertices to reconstruct a 3D vertex" << std::endl;

        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            if (canvases[view])
                delete canvases[view];
        }

        return STATUS_CODE_NOT_FOUND;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::MakeNetworkInputFromHits(const CaloHitList &caloHits, const HitType view, const float xMin,
    const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float driftStep{0.5f};

    // Determine the bin edges
    std::vector<double> xBinEdges(m_width + 1);
    std::vector<double> zBinEdges(m_height + 1);
    xBinEdges[0] = xMin - 0.5f * driftStep;
    const double dx = ((xMax + 0.5f * driftStep) - xBinEdges[0]) / m_width;
    for (int i = 1; i < m_width + 1; ++i)
        xBinEdges[i] = xBinEdges[i - 1] + dx;
    zBinEdges[0] = zMin - 0.5f * pitch;
    const double dz = ((zMax + 0.5f * pitch) - zBinEdges[0]) / m_height;
    for (int i = 1; i < m_height + 1; ++i)
        zBinEdges[i] = zBinEdges[i - 1] + dz;

    LArDLHelper::InitialiseInput({1, 1, m_height, m_width}, networkInput);
    auto accessor = networkInput.accessor<float, 4>();

    for (const CaloHit *pCaloHit : caloHits)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        if (m_pass > 1)
        {
            if (x < xMin || x > xMax || z < zMin || z > zMax)
                continue;
        }
        const float adc{pCaloHit->GetMipEquivalentEnergy()};
        const int pixelX{static_cast<int>(std::floor((x - xBinEdges[0]) / dx))};
        const int pixelZ{static_cast<int>(std::floor((z - zBinEdges[0]) / dz))};
        accessor[0][0][pixelZ][pixelX] += adc;
    }
    for (int row = 0; row < m_height; ++row)
    {
        for (int col = 0; col < m_width; ++col)
        {
            const float value{accessor[0][0][row][col]};
            if (value > 0)
                pixelVector.emplace_back(std::make_pair(row, col));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::GetNetworkVertices(const CanvasViewMap &canvases, CartesianPointVector &positionVector) const
{
    CartesianPointVector verticesU, verticesV, verticesW;
    this->GetVerticesFromCanvas(*canvases.at(TPC_VIEW_U), verticesU);
    this->GetVerticesFromCanvas(*canvases.at(TPC_VIEW_V), verticesV);
    this->GetVerticesFromCanvas(*canvases.at(TPC_VIEW_W), verticesW);

    std::vector<VertexTuple> vertexTuples;
    int nEmptyLists{(verticesU.empty() ? 1 : 0) + (verticesV.empty() ? 1 : 0) + (verticesW.empty() ? 1 : 0)};
    
    if (nEmptyLists == 0)
    {
        for (const CartesianVector &vertexU : verticesU)
        {
            for (const CartesianVector &vertexV : verticesV)
            {
                for (const CartesianVector &vertexW : verticesW)
                {
                    const float xMin{std::min({vertexU.GetX(), vertexV.GetX(), vertexW.GetX()})};
                    const float xMax{std::max({vertexU.GetX(), vertexV.GetX(), vertexW.GetX()})};
                    if ((xMax - xMin) < 10.f)
                    {
                        const VertexTuple &tuple{VertexTuple(this->GetPandora(), vertexU, vertexV, vertexW)};
                        if (tuple.GetChi2() < 1)
                            vertexTuples.emplace_back(tuple);
                    }
                }
            }
        }
    }
    else if (nEmptyLists == 1)
    {
        if (verticesU.empty())
        {
            for (const CartesianVector &vertexV : verticesV)
            {
                for (const CartesianVector &vertexW : verticesW)
                {
                    const float xMin{std::min({vertexV.GetX(), vertexW.GetX()})};
                    const float xMax{std::max({vertexV.GetX(), vertexW.GetX()})};
                    if ((xMax - xMin) < 10.f)
                    {
                        const VertexTuple &tuple{VertexTuple(this->GetPandora(), vertexV, vertexW, TPC_VIEW_V, TPC_VIEW_W)};
                        if (tuple.GetChi2() < 1)
                            vertexTuples.emplace_back(tuple);
                    }
                }
            }
        }
        else if (verticesV.empty())
        {
            for (const CartesianVector &vertexU : verticesU)
            {
                for (const CartesianVector &vertexW : verticesW)
                {
                    const float xMin{std::min({vertexU.GetX(), vertexW.GetX()})};
                    const float xMax{std::max({vertexU.GetX(), vertexW.GetX()})};
                    if ((xMax - xMin) < 10.f)
                    {
                        const VertexTuple &tuple{VertexTuple(this->GetPandora(), vertexU, vertexW, TPC_VIEW_U, TPC_VIEW_W)};
                        if (tuple.GetChi2() < 1)
                            vertexTuples.emplace_back(tuple);
                    }
                }
            }
        }
        else
        {
            for (const CartesianVector &vertexU : verticesU)
            {
                for (const CartesianVector &vertexV : verticesV)
                {
                    const float xMin{std::min({vertexU.GetX(), vertexV.GetX()})};
                    const float xMax{std::max({vertexU.GetX(), vertexV.GetX()})};
                    if ((xMax - xMin) < 10.f)
                    {
                        const VertexTuple &tuple{VertexTuple(this->GetPandora(), vertexU, vertexV, TPC_VIEW_U, TPC_VIEW_V)};
                        if (tuple.GetChi2() < 1)
                            vertexTuples.emplace_back(tuple);
                    }
                }
            }
        }
    }
    else
    {
        std::cout << "Insufficient 2D vertices to reconstruct a 3D vertex" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    // Sort the vertex tuples here and pick the unique ones that look sound
    std::sort(vertexTuples.begin(), vertexTuples.end(),
        [](const VertexTuple &tuple1, const VertexTuple &tuple2)
        {
            const CartesianPointVector &components1{tuple1.GetComponents()};
            const CartesianPointVector &components2{tuple2.GetComponents()};
            if (components1.size() == components2.size())
                return tuple1.GetChi2() < tuple2.GetChi2();
            else
                return components1.size() > components2.size();
        });

    CartesianPointVector used;
    for (const VertexTuple &tuple : vertexTuples)
    {
        const CartesianPointVector &components{tuple.GetComponents()};
        bool isAvailable{true};
        for (const CartesianVector &component : components)
        {
            if (std::find(used.begin(), used.end(), component) != used.end())
            {
                isAvailable = false;
                break;
            }
        }
        if (!isAvailable || (components.size() < 3 && tuple.GetChi2() > 0.1))
            continue;

        positionVector.emplace_back(tuple.GetPosition());
        for (const CartesianVector &component : components)
            used.emplace_back(component);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::GetVerticesFromCanvas(const Canvas &canvas, CartesianPointVector &vertices) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);

    HitType view{canvas.m_view};
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float driftStep{0.5f};

    const double dx{((canvas.m_xMax + 0.5f * driftStep) - (canvas.m_xMin - 0.5f * driftStep)) / m_width};
    const double dz{((canvas.m_zMax + 0.5f * pitch) - (canvas.m_zMin - 0.5f * pitch)) / m_height};

    float maxIntensity{0.f};
    for (int xp = 0; xp < m_width; ++xp)
    {
        int xpp{xp + canvas.m_colOffset};
        if (xpp >= canvas.m_width)
            continue;

        for (int zp = 0; zp < m_height; ++zp)
        {
            int zpp{zp + canvas.m_rowOffset};
            if (zpp >= canvas.m_height)
                continue;

            const float localIntensity{canvas.m_canvas[zpp][xpp]};
            if (localIntensity > maxIntensity)
                maxIntensity = localIntensity;
        }
    }
    const float threshold{maxIntensity * 0.3f};

    std::vector<std::vector<std::pair<int, int>>> peaks;
    for (int xp = 0; xp < m_width; ++xp)
    {
        int xpp{xp + canvas.m_colOffset};
        if (xpp >= canvas.m_width)
            continue;

        for (int zp = 0; zp < m_height; ++zp)
        {
            int zpp{zp + canvas.m_rowOffset};
            if (zpp >= canvas.m_height)
                continue;
            const float localIntensity{canvas.m_canvas[zpp][xpp]};

            std::vector<std::pair<int, int>> peak;
            bool hasLowNeighbour{false};
            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <=1; ++dc)
                {
                    if (dr == 0 && dc == 0)
                        continue;
                    const int r{zpp + dr}, c{xpp + dc};
                    if (r < 0 || r >= canvas.m_height || c < 0 || c >= canvas.m_width)
                        continue;

                    const float neighborIntensity{canvas.m_canvas[r][c]};
                    if (localIntensity > neighborIntensity)
                    {
                        hasLowNeighbour = true;
                    }
                    else if (localIntensity < neighborIntensity)
                    {
                        hasLowNeighbour = false;
                        break;
                    }
                }
            }
            if (hasLowNeighbour && localIntensity > threshold)
                this->GrowPeak(canvas, xpp, zpp, localIntensity, peak);
            if (!peak.empty())
                peaks.emplace_back(peak);
        }
    }

    for (const auto &peak : peaks)
    {
        float row{0}, col{0};
        for (const auto &pixel : peak)
        {
            row += pixel.second - canvas.m_rowOffset;
            col += pixel.first - canvas.m_colOffset;
        }
        row /= peak.size();
        col /= peak.size();

        const float x{static_cast<float>(col * dx + canvas.m_xMin)};
        const float z{static_cast<float>(row * dz + canvas.m_zMin)};
        CartesianVector pt(x, 0, z);
        vertices.emplace_back(pt);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DlSecondaryVertexingAlgorithm::GrowPeak(const Canvas &canvas, int col, int row, float intensity, std::vector<std::pair<int, int>> &peak) const
{
//    if (col >= 0 && col < canvas.m_width && row >= 0 && row < canvas.m_height)
//    {
//        std::cout << "Visited: " << canvas.m_visited[row][col] << " " << 
//    }
    if (col < 0 || col >= canvas.m_width || row < 0 || row >= canvas.m_height || canvas.m_visited[row][col] || canvas.m_canvas[row][col] < intensity)
        return false;

    // Check that no adjacent pixel is larger than this one
    for (int dc = -1; dc <= 1; ++dc)
    {
        const int c{col + dc};
        if (c < 0 || c >= canvas.m_width)
            continue;

        for (int dr = -1; dr <= 1; ++dr)
        {
            const int r{row + dr};
            if (r < 0 || r >= canvas.m_height)
                continue;

            if (dr == 0 && dc == 0)
                continue;

            const float neighborIntensity{canvas.m_canvas[r][c]};
            if (neighborIntensity > intensity)
                return false;
        }
    }

    // Need to check we aren't growing into a higher peak, if we are restart from the current pixel
    float localIntensity{canvas.m_canvas[row][col]};
    if (localIntensity > intensity)
    {
        intensity = localIntensity;
        for (const auto &pixel : peak)
            canvas.m_visited[pixel.second][pixel.first] = false;
        peak.clear();
        //std::cout << ". New size " << peak.size() << " new intensity " << intensity << std::endl;
        this->GrowPeak(canvas, col, row, intensity, peak);
        return true;
    }

    // Add pixel to the peak
    canvas.m_visited[row][col] = true;
    peak.emplace_back(std::make_pair(col, row));
    //std::cout << "   Added" << std::endl;

    for (int dc = -1; dc <= 1; ++dc)
    {
        for (int dr = -1; dr <= 1; ++dr)
        {
            if (dr == 0 && dc == 0)
                continue;
            //std::cout << "   Adjacent (" << (row + i) << "," << (col + j) << ")" << std::endl;
            bool reset{this->GrowPeak(canvas, col + dc, row + dr, intensity, peak)};
            // If we started growing a non-peak region, stop looking relative to the previous peak
            if (reset)
                return reset;
        }
    }

    return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSecondaryVertexingAlgorithm::GetCanvasParameters(const LArDLHelper::TorchOutput &networkOutput, const PixelVector &pixelVector,
    int &colOffset, int &rowOffset, int &width, int &height) const
{
    const double scaleFactor{std::sqrt(m_height * m_height + m_width * m_width)};
    // output is a 1 x num_classes x height x width tensor
    // we want the maximum value in the num_classes dimension (1) for every pixel
    auto classes{torch::argmax(networkOutput, 1)};
    // the argmax result is a 1 x height x width tensor where each element is a class id
    auto classesAccessor{classes.accessor<long, 3>()};
    int colOffsetMin{0}, colOffsetMax{0}, rowOffsetMin{0}, rowOffsetMax{0};
    for (const auto &[row, col] : pixelVector)
    {
        const auto cls{classesAccessor[0][row][col]};
        const double threshold{m_thresholds[cls]};
        if (threshold > 0. && threshold < 1.)
        {
            const int distance = static_cast<int>(std::round(std::ceil(scaleFactor * threshold)));
            if ((row - distance) < rowOffsetMin)
                rowOffsetMin = row - distance;
            if ((row + distance) > rowOffsetMax)
                rowOffsetMax = row + distance;
            if ((col - distance) < colOffsetMin)
                colOffsetMin = col - distance;
            if ((col + distance) > colOffsetMax)
                colOffsetMax = col + distance;
        }
    }
    colOffset = colOffsetMin < 0 ? -colOffsetMin : 0;
    rowOffset = rowOffsetMin < 0 ? -rowOffsetMin : 0;
    width = std::max(colOffsetMax + colOffset + 1, m_width);
    height = std::max(rowOffsetMax + rowOffset + 1, m_height);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSecondaryVertexingAlgorithm::DrawRing(float **canvas, const int row, const int col, const int inner, const int outer, const float weight) const
{
    // Set the starting position for each circle bounding the ring
    int c1{inner}, r1{0}, c2{outer}, r2{0};
    int inner2{inner * inner}, outer2{outer * outer};
    while (c2 >= r2)
    {
        // Set the output pixel location
        int rp2{r2}, cp2{c2};
        // We're still within the octant for the inner ring, so use the inner pixel location (see Update comment below)
        // Note also that the inner row is always the same as the outer row, so no need to define rp1
        int cp1{c1};
        if (c1 <= r1)
        { // We've completed the arc of the inner ring already, so just move radially out from here (see Update comment below)
            cp1 = r2;
        }
        // Fill the pixels from inner to outer in the current row and their mirror pixels in the other octants
        for (int c = cp1; c <= cp2; ++c)
        {
            canvas[row + rp2][col + c] += weight;
            if (rp2 != c)
                canvas[row + c][col + rp2] += weight;
            if (rp2 != 0 && cp2 != 0)
            {
                canvas[row - rp2][col - c] += weight;
                if (rp2 != c)
                    canvas[row - c][col - rp2] += weight;
            }
            if (rp2 != 0)
            {
                canvas[row - rp2][col + c] += weight;
                if (rp2 != c)
                    canvas[row + c][col - rp2] += weight;
            }
            if (cp2 != 0)
            {
                canvas[row + rp2][col - c] += weight;
                if (rp2 != c)
                    canvas[row - c][col + rp2] += weight;
            }
        }
        // Only update the inner location while it remains in the octant (outer ring also remains in the octant of course, but the logic of
        // the update means that the inner ring can leave its octant before the outer ring is complete, so we need to stop that)
        if (c1 > r1)
            this->Update(inner2, c1, r1);
        // Update the outer location - increase the row position with every step, decrease the column position if conditions are met
        this->Update(outer2, c2, r2);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSecondaryVertexingAlgorithm::Update(const int radius2, int &col, int &row) const
{
    // Bresenham midpoint circle algorithm to determine if we should update the column position
    // This obscure looking block of code uses bit shifts and integer arithmetic to perform this check as efficiently as possible
    const int a{1 - (col << 2)};
    const int b{col * col + row * row - radius2 + (row << 2) + 1};
    const int c{(a << 2) * b + a * a};
    if (c < 0)
    {
        --col;
        ++row;
    }
    else
        ++row;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList2D, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCParticleList &mcHierarchy) const
{
    try
    {
        for (const auto &[mc, hits] : mcToHitsMap)
        {
            (void)hits;
            mcHierarchy.push_back(mc);
            LArMCParticleHelper::GetAllAncestorMCParticles(mc, mcHierarchy);
        }
    }
    catch (const StatusCodeException &e)
    {
        return e.GetStatusCode();
    }

    // Move the neutrino to the front of the list
    auto pivot =
        std::find_if(mcHierarchy.begin(), mcHierarchy.end(), [](const MCParticle *mc) -> bool { return LArMCParticleHelper::IsNeutrino(mc); });
    (void)pivot;
    if (pivot != mcHierarchy.end())
        std::rotate(mcHierarchy.begin(), pivot, std::next(pivot));
    else
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSecondaryVertexingAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
    // Find the range of x and z values in the view
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }

    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType view{caloHitList.front()->GetHitType()};
    const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
    if (!(isU || isV || isW))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

    // Avoid unreasonable rescaling of very small hit regions, pixels are assumed to be 0.5cm in x and wire pitch in z
    // ATTN: Rescaling is to a size 1 pixel smaller than the intended image to ensure all hits fit within an imaged binned
    // to be one pixel wider than this
    const float xRange{xMax - xMin}, zRange{zMax - zMin};
    const float minXSpan{m_driftStep * (m_width - 1)};
    if (xRange < minXSpan)
    {
        const float padding{0.5f * (minXSpan - xRange)};
        xMin -= padding;
        xMax += padding;
    }
    const float minZSpan{pitch * (m_height - 1)};
    if (zRange < minZSpan)
    {
        const float padding{0.5f * (minZSpan - zRange)};
        zMin -= padding;
        zMax += padding;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::MakeCandidateVertexList(const CartesianPointVector &positions)
{
    const VertexList *pVertexList{nullptr};
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const CartesianVector &position : positions)
    {
        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pVertex(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSecondaryVertexingAlgorithm::GetTrueVertexPosition(float &x, float &y, float &z) const
{
    const CartesianVector &trueVertex{this->GetTrueVertex()};
    x = trueVertex.GetX();
    y = trueVertex.GetY();
    z = trueVertex.GetZ();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSecondaryVertexingAlgorithm::GetTrueVertexPosition(float &x, float &u, float &v, float &w) const
{
    const CartesianVector &trueVertex{this->GetTrueVertex()};
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    x = trueVertex.GetX();
    u = static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ()));
    v = static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ()));
    w = static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ()));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &DlSecondaryVertexingAlgorithm::GetTrueVertex() const
{
    const MCParticleList *pMCParticleList{nullptr};
    if (STATUS_CODE_SUCCESS == PandoraContentApi::GetCurrentList(*this, pMCParticleList) && pMCParticleList)
    {
        MCParticleVector primaries;
        LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
        if (!primaries.empty())
        {
            const MCParticle *primary{primaries.front()};
            const MCParticleList &parents{primary->GetParentList()};
            if (parents.size() == 1)
            {
                const MCParticle *trueNeutrino{parents.front()};
                if (LArMCParticleHelper::IsNeutrino(trueNeutrino))
                    return primaries.front()->GetVertex();
            }
        }
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSecondaryVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Pass", m_pass));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_height));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_width));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DistanceThresholds", m_thresholds));
    m_nClasses = m_thresholds.size() - 1;
    if (m_pass > 1)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));
    }

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
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
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
        if (m_writeTree)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        }
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

DlSecondaryVertexingAlgorithm::Canvas::Canvas(const HitType view, const int width, const int height, const int colOffset, const int rowOffset,
    const float xMin, const float xMax, const float zMin, const float zMax) :
    m_view{view},
    m_width{width},
    m_height{height},
    m_colOffset{colOffset},
    m_rowOffset{rowOffset},
    m_xMin{xMin},
    m_xMax{xMax},
    m_zMin{zMin},
    m_zMax{zMax}
{
    m_canvas = new float*[m_height];
    m_visited = new bool*[m_height];
    for (int r = 0; r < m_height; ++r)
    {
        m_canvas[r] = new float[m_width]{};
        m_visited[r] = new bool[m_width]{false};
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlSecondaryVertexingAlgorithm::Canvas::~Canvas()
{
    for (int r = 0; r < m_height; ++r)
    {
        delete[] m_canvas[r];
        delete[] m_visited[r];
    }
    delete[] m_canvas;
    delete[] m_visited;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

DlSecondaryVertexingAlgorithm::VertexTuple::VertexTuple(const Pandora &pandora, const CartesianVector &vertexU, const CartesianVector &vertexV,
    const CartesianVector &vertexW) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    LArGeometryHelper::MergeThreePositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vertexU, vertexV, vertexW, m_pos, m_chi2);
    if (m_chi2 <= 1.f)
    {
        m_components.emplace_back(vertexU);
        m_components.emplace_back(vertexV);
        m_components.emplace_back(vertexW);
    }
    else
    {
        CartesianVector vertexUV(0.f, 0.f, 0.f);
        float chi2UV{0.f};
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, vertexU, vertexV, vertexUV, chi2UV);

        CartesianVector vertexUW(0.f, 0.f, 0.f);
        float chi2UW{0.f};
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_U, TPC_VIEW_W, vertexU, vertexW, vertexUW, chi2UW);

        CartesianVector vertexVW(0.f, 0.f, 0.f);
        float chi2VW{0.f};
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_V, TPC_VIEW_W, vertexV, vertexW, vertexVW, chi2VW);
         
        if (chi2UV < m_chi2)
        {
            m_pos = vertexUV;
            m_chi2 = chi2UV;
            m_components.emplace_back(vertexU);
            m_components.emplace_back(vertexV);
        }
        if (chi2UW < m_chi2)
        {
            m_pos = vertexUW;
            m_chi2 = chi2UW;
            m_components.clear();
            m_components.emplace_back(vertexU);
            m_components.emplace_back(vertexW);
        }
        if (chi2VW < m_chi2)
        {
            m_pos = vertexVW;
            m_chi2 = chi2VW;
            m_components.clear();
            m_components.emplace_back(vertexV);
            m_components.emplace_back(vertexW);
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlSecondaryVertexingAlgorithm::VertexTuple::VertexTuple(const Pandora &pandora, const CartesianVector &vertex1, const CartesianVector &vertex2,
    const HitType view1, const HitType view2) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    m_components.emplace_back(vertex1);
    m_components.emplace_back(vertex2);

    LArGeometryHelper::MergeTwoPositions3D(pandora, view1, view2, vertex1, vertex2, m_pos, m_chi2);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &DlSecondaryVertexingAlgorithm::VertexTuple::GetPosition() const
{
    return m_pos;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianPointVector &DlSecondaryVertexingAlgorithm::VertexTuple::GetComponents() const
{
    return m_components;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float DlSecondaryVertexingAlgorithm::VertexTuple::GetChi2() const
{
    return m_chi2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::string DlSecondaryVertexingAlgorithm::VertexTuple::ToString() const
{
    const float x{m_pos.GetX()}, y{m_pos.GetY()}, z{m_pos.GetZ()};
    return "3D pos: (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")   X2 = " + std::to_string(m_chi2);
}

} // namespace lar_dl_content