/**
 *  @file   larpandoradlcontent/LArVertex/DlVertexingAlgorithm.cc
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

#include "larpandoradlcontent/LArVertex/DlVertexingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlVertexingAlgorithm::DlVertexingAlgorithm() :
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

DlVertexingAlgorithm::~DlVertexingAlgorithm()
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

StatusCode DlVertexingAlgorithm::Run()
{
    ++m_event;

    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlVertexingAlgorithm::PrepareTrainingSample()
{
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    LArEventTopology eventTopology(*pCaloHitList2D);
    eventTopology.ConstructVisibleHierarchy();
    eventTopology.PruneHierarchy();
    CartesianPointVector vertices;
    eventTopology.GetVertices(vertices);

    if (vertices.empty())
        return STATUS_CODE_SUCCESS;

/*    for (const CartesianVector &vertex : vertices)
    {
        std::string str{"(" + std::to_string(vertex.GetX()) + "," + std::to_string(vertex.GetY()) + "," + std::to_string(vertex.GetZ()) + ")"};
        std::cout << "(" << vertex.GetX() << "," << vertex.GetY() << "," << vertex.GetZ() << ")" << std::endl;
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertex, str, BLUE, 2));
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMCToHitsMap(mcToHitsMap));
    MCParticleList hierarchy;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CompleteMCHierarchy(mcToHitsMap, hierarchy));

    // Get boundaries for hits and make x dimension common

    float xMin{0.f}, xMax{0.f}, yMin{0.f}, yMax{0.f}, zMin{0.f}, zMax{0.f};
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList3D", pCaloHitList));
    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    this->GetHitRegion(*pCaloHitList, xMin, xMax, yMin, yMax, zMin, zMax);

    const std::string trainingFilename{m_trainingOutputFile + ".csv"};
    const unsigned long nVertices{vertices.size()};
    unsigned long nHits{0};

    // Only train on events where there is a vertex in the fiducial volume
    bool hasFiducialVertex{false};
    for (const CartesianVector vertex : vertices)
    {
        if (LArVertexHelper::IsInFiducialVolume(this->GetPandora(), vertex, "dune_fd_hd"))
        {
            hasFiducialVertex = true;
            break;
        }
    }

    if (!hasFiducialVertex)
        return STATUS_CODE_SUCCESS;

    LArMvaHelper::MvaFeatureVector featureVector;
    featureVector.emplace_back(static_cast<double>(m_event));
    featureVector.emplace_back(static_cast<double>(nVertices));
    // Vertices
    for (const CartesianVector &vertex : vertices)
    {
        featureVector.emplace_back(vertex.GetX());
        featureVector.emplace_back(vertex.GetY());
        featureVector.emplace_back(vertex.GetZ());
    }

    // Retain the hit region
    featureVector.emplace_back(xMin);
    featureVector.emplace_back(xMax);
    featureVector.emplace_back(yMin);
    featureVector.emplace_back(yMax);
    featureVector.emplace_back(zMin);
    featureVector.emplace_back(zMax);

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        // If on a refinement pass, drop hits outside the region of interest
        // NEEDS RECONSIDERATION IN THE CONTEXT OF ORTHOGONAL VIEWS
        //if (m_pass > 1 && (x < xMin || x > xMax || z < zMin || z > zMax))
        //    continue;
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetY()));
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetMipEquivalentEnergy()));
        ++nHits;
    }
    featureVector.insert(featureVector.begin() + 2 + 3 * nVertices + 6, static_cast<double>(nHits));
    // Only write out the feature vector if there were enough hits in the region of interest
    if (nHits > 10)
        LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::Infer()
{
    float xMin{0.f}, xMax{0.f}, yMin{0.f}, yMax{0.f}, zMin{0.f}, zMax{0.f};
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList3D", pCaloHitList));
    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    this->GetHitRegion(*pCaloHitList, xMin, xMax, yMin, yMax, zMin, zMax);

    CanvasViewMap canvases;
    CartesianPointVector vertexCandidatesU, vertexCandidatesV, vertexCandidatesW;
    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        LArDLHelper::TorchInput input;
        PixelVector pixelVector;

        switch (view)
        {
            case TPC_VIEW_U:
                this->MakeNetworkInputFromHits(*pCaloHitList, view, xMin, xMax, yMin, yMax, input, pixelVector);
                break;
            case TPC_VIEW_V:
                this->MakeNetworkInputFromHits(*pCaloHitList, view, xMin, xMax, zMin, zMax, input, pixelVector);
                break;
            default:
                this->MakeNetworkInputFromHits(*pCaloHitList, view, yMin, yMax, zMin, zMax, input, pixelVector);
                break;
        }

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
        std::cout << "View " << view << " " << canvasWidth << " " << canvasHeight << " " << colOffset << " " << rowOffset << std::endl;

        switch (view)
        {
            case TPC_VIEW_U:
                canvases[view] = new Canvas(canvasWidth, canvasHeight, colOffset, rowOffset, xMin, xMax, yMin, yMax);
                break;
            case TPC_VIEW_V:
                canvases[view] = new Canvas(canvasWidth, canvasHeight, colOffset, rowOffset, xMin, xMax, zMin, zMax);
                break;
            default:
                canvases[view] = new Canvas(canvasWidth, canvasHeight, colOffset, rowOffset, yMin, yMax, zMin, zMax);
                break;
        }

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

        LArMvaHelper::MvaFeatureVector featureVector;
        std::cout << view << " " << canvases[view]->m_width << " " << canvases[view]->m_height << " " << canvases[view]->m_colOffset << " " <<
            canvases[view]->m_rowOffset << std::endl;
        featureVector.emplace_back(static_cast<double>(canvases[view]->m_width));
        featureVector.emplace_back(static_cast<double>(canvases[view]->m_height));
        for (const Pixel &pixel : pixelVector)
        {
            featureVector.emplace_back(static_cast<double>(pixel.first + canvases[view]->m_rowOffset));
            featureVector.emplace_back(static_cast<double>(pixel.second + canvases[view]->m_colOffset));
        }

        LArMvaHelper::ProduceTrainingExample("hits.csv", true, featureVector);
    }

    for (HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        LArMvaHelper::MvaFeatureVector featureVector;
        std::cout << view << " " << canvases[view]->m_width << " " << canvases[view]->m_height << " " << canvases[view]->m_colOffset << " " <<
            canvases[view]->m_rowOffset << std::endl;
        featureVector.emplace_back(static_cast<double>(canvases[view]->m_width));
        featureVector.emplace_back(static_cast<double>(canvases[view]->m_height));
        featureVector.emplace_back(static_cast<double>(canvases[view]->m_colOffset));
        featureVector.emplace_back(static_cast<double>(canvases[view]->m_rowOffset));
        for (int row = 0; row < canvases[view]->m_height; ++row)
            for (int col = 0; col < canvases[view]->m_width; ++col)
                featureVector.emplace_back(static_cast<double>(canvases[view]->m_canvas[row][col]));

        LArMvaHelper::ProduceTrainingExample("canvases.csv", true, featureVector);
    }

    CartesianPointVector vertexVector;
    this->MakeWirePlaneCoordinatesFromCanvas(canvases, vertexVector);

    if (!vertexVector.empty())
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

        for (const CartesianVector &vertex : vertexVector)
        {
            const CartesianVector vec2d(vertex.GetX(), 0, vertex.GetZ());
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vec2d, "vec", BLUE, 1));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        std::cout << "Final: " << vertexVector.front() << std::endl;
        StatusCode status{this->MakeCandidateVertexList(vertexVector)};

        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            delete canvases[view];

        return status;
    }
    else
    {
        std::cout << "Insufficient 2D vertices to reconstruct a 3D vertex" << std::endl;

        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            delete canvases[view];

        return STATUS_CODE_NOT_FOUND;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeNetworkInputFromHits(const CaloHitList &caloHits, const HitType view, const float xMin,
    const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const
{
    const float pitch(0.5f);
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
        const float x{view == TPC_VIEW_W ? pCaloHit->GetPositionVector().GetY() : pCaloHit->GetPositionVector().GetX()};
        const float z{view == TPC_VIEW_U ? pCaloHit->GetPositionVector().GetY() : pCaloHit->GetPositionVector().GetZ()};
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

StatusCode DlVertexingAlgorithm::MakeWirePlaneCoordinatesFromCanvas(const CanvasViewMap &canvases, CartesianPointVector &positionVector) const
{
    const float pitch(0.5f);
    const float driftStep{0.5f};

    const double dx = ((canvases.at(TPC_VIEW_U)->m_xMax + 0.5f * driftStep) - (canvases.at(TPC_VIEW_U)->m_xMin - 0.5f * driftStep)) / m_width;
    const double dy = ((canvases.at(TPC_VIEW_U)->m_zMax + 0.5f * driftStep) - (canvases.at(TPC_VIEW_U)->m_zMin - 0.5f * driftStep)) / m_width;
    const double dz = ((canvases.at(TPC_VIEW_V)->m_zMax + 0.5f * pitch) - (canvases.at(TPC_VIEW_V)->m_zMin - 0.5f * pitch)) / m_height;

    float best{-1.f};
    int xBest{0}, yBest{0}, zBest{0};
/*    for (int xp = 0; xp < m_width; ++xp)
    {
        if ((xp + canvases.at(TPC_VIEW_U)->m_colOffset >= canvases.at(TPC_VIEW_U)->m_width) ||
            (xp + canvases.at(TPC_VIEW_V)->m_colOffset >= canvases.at(TPC_VIEW_V)->m_width))
            continue;

        for (int yp = 0; yp < m_width; ++yp)
        {
            if ((yp + canvases.at(TPC_VIEW_U)->m_rowOffset >= canvases.at(TPC_VIEW_U)->m_height) ||
                (yp + canvases.at(TPC_VIEW_W)->m_colOffset >= canvases.at(TPC_VIEW_W)->m_width))
                continue;

            for (int zp = 0; zp < m_width; ++zp)
            {
                if ((zp + canvases.at(TPC_VIEW_V)->m_rowOffset >= canvases.at(TPC_VIEW_V)->m_height) ||
                    (zp + canvases.at(TPC_VIEW_W)->m_rowOffset >= canvases.at(TPC_VIEW_W)->m_height))
                    continue;

                float contribution{canvases.at(TPC_VIEW_U)->m_canvas[yp + canvases.at(TPC_VIEW_U)->m_rowOffset][xp + canvases.at(TPC_VIEW_U)->m_colOffset]};
                contribution += canvases.at(TPC_VIEW_V)->m_canvas[yp + canvases.at(TPC_VIEW_V)->m_rowOffset][xp + canvases.at(TPC_VIEW_V)->m_colOffset];
                contribution += canvases.at(TPC_VIEW_W)->m_canvas[zp + canvases.at(TPC_VIEW_W)->m_rowOffset][yp + canvases.at(TPC_VIEW_W)->m_colOffset];
                if (contribution > best)
                {
                    best = contribution;
                    xBest = xp;
                    yBest = yp;
                    zBest = zp;
                }
            }
        }
    }*/

    for (int xp = 0; xp < m_width; ++xp)
    {
        if (xp + canvases.at(TPC_VIEW_V)->m_colOffset >= canvases.at(TPC_VIEW_V)->m_width)
            continue;

        for (int zp = 0; zp < m_width; ++zp)
        {
            if (zp + canvases.at(TPC_VIEW_V)->m_rowOffset >= canvases.at(TPC_VIEW_V)->m_height)
                continue;

            float contribution{canvases.at(TPC_VIEW_V)->m_canvas[zp + canvases.at(TPC_VIEW_V)->m_rowOffset][xp + canvases.at(TPC_VIEW_V)->m_colOffset]};
            if (contribution > best)
            {
                best = contribution;
                xBest = xp;
                yBest = 0;
                zBest = zp;
            }
        }
    }


    const float x{static_cast<float>(xBest * dx + canvases.at(TPC_VIEW_U)->m_xMin)};
    const float y{static_cast<float>(yBest * dy + canvases.at(TPC_VIEW_U)->m_zMin)};
    const float z{static_cast<float>(zBest * dz + canvases.at(TPC_VIEW_V)->m_zMin)};

    CartesianVector pt(x, y, z);
    positionVector.emplace_back(pt);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::GetCanvasParameters(const LArDLHelper::TorchOutput &networkOutput, const PixelVector &pixelVector,
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

void DlVertexingAlgorithm::DrawRing(float **canvas, const int row, const int col, const int inner, const int outer, const float weight) const
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

void DlVertexingAlgorithm::Update(const int radius2, int &col, int &row) const
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

StatusCode DlVertexingAlgorithm::GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
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

StatusCode DlVertexingAlgorithm::CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCParticleList &mcHierarchy) const
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

void DlVertexingAlgorithm::GetHitRegion(const CaloHitList &caloHitList, HitType view, float &xMin, float &xMax, float &zMin, float &zMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
    // Find the range of x and z values in the view
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{view == TPC_VIEW_W ? pCaloHit->GetPositionVector().GetY() : pCaloHit->GetPositionVector().GetX()};
        const float z{view == TPC_VIEW_U ? pCaloHit->GetPositionVector().GetY() : pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }

    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
    if (!(isU || isV || isW))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    const float pitch(0.5f);

    if (m_pass > 1)
    {
        const VertexList *pVertexList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pVertexList));
        if (pVertexList->empty())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const CartesianVector &vertex{pVertexList->front()->GetPosition()};

        // Get hit distribution left/right asymmetry
        int nHitsLeft{0}, nHitsRight{0};
        const double xVtx{vertex.GetX()};
        for (const std::string &listname : m_caloHitListNames)
        {
            const CaloHitList *pCaloHitList(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
            if (pCaloHitList->empty())
                continue;
            for (const CaloHit *const pCaloHit : *pCaloHitList)
            {
                const CartesianVector &pos{pCaloHit->GetPositionVector()};
                if (pos.GetX() <= xVtx)
                    ++nHitsLeft;
                else
                    ++nHitsRight;
            }
        }
        const int nHitsTotal{nHitsLeft + nHitsRight};
        if (nHitsTotal == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const float xAsymmetry{nHitsLeft / static_cast<float>(nHitsTotal)};

        // Vertices
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        double zVtx{0.};
        if (isW)
            zVtx += transform->YZtoW(vertex.GetY(), vertex.GetZ());
        else if (isV)
            zVtx += transform->YZtoV(vertex.GetY(), vertex.GetZ());
        else
            zVtx = transform->YZtoU(vertex.GetY(), vertex.GetZ());

        // Get hit distribution upstream/downstream asymmetry
        int nHitsUpstream{0}, nHitsDownstream{0};
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            if (pos.GetZ() <= zVtx)
                ++nHitsUpstream;
            else
                ++nHitsDownstream;
        }
        const int nHitsViewTotal{nHitsUpstream + nHitsDownstream};
        if (nHitsViewTotal == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const float zAsymmetry{nHitsUpstream / static_cast<float>(nHitsViewTotal)};

        const float xSpan{m_driftStep * (m_width - 1)};
        xMin = xVtx - xAsymmetry * xSpan;
        xMax = xMin + (m_driftStep * (m_width - 1));
        const float zSpan{pitch * (m_height - 1)};
        zMin = zVtx - zAsymmetry * zSpan;
        zMax = zMin + zSpan;
    }

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

void DlVertexingAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &yMin, float &yMax, float &zMin,
    float &zMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    yMin = std::numeric_limits<float>::max();
    yMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();

    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Find the range of x and z values in the view
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float y{pCaloHit->GetPositionVector().GetY()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        yMin = std::min(y, yMin);
        yMax = std::max(y, yMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }

    const float pitch(0.5f);

    if (m_pass > 1)
    {
        // NOT YET IMPLEMENTED FOR ORTHOGONAL VIEWS
    }

    // Avoid unreasonable rescaling of very small hit regions, pixels are assumed to be 0.5cm in x and wire pitch in z
    // ATTN: Rescaling is to a size 1 pixel smaller than the intended image to ensure all hits fit within an imaged binned
    // to be one pixel wider than this
    const float xRange{xMax - xMin};
    const float minXSpan{m_driftStep * (m_width - 1)};
    if (xRange < minXSpan)
    {
        const float padding{0.5f * (minXSpan - xRange)};
        xMin -= padding;
        xMax += padding;
    }
    const float yRange{yMax - yMin};
    const float minYSpan{pitch * (m_width - 1)};
    if (yRange < minYSpan)
    {
        const float padding{0.5f * (minYSpan - yRange)};
        yMin -= padding;
        yMax += padding;
    }
    const float zRange{zMax - zMin};
    const float minZSpan{pitch * (m_height - 1)};
    if (zRange < minZSpan)
    {
        const float padding{0.5f * (minZSpan - zRange)};
        zMin -= padding;
        zMax += padding;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::MakeCandidateVertexList(const CartesianPointVector &positions)
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
        std::cout << "All good" << std::endl;
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::GetTrueVertexPosition(float &x, float &y, float &z) const
{
    const CartesianVector &trueVertex{this->GetTrueVertex()};
    x = trueVertex.GetX();
    y = trueVertex.GetY();
    z = trueVertex.GetZ();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::GetTrueVertexPosition(float &x, float &u, float &v, float &w) const
{
    const CartesianVector &trueVertex{this->GetTrueVertex()};
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    x = trueVertex.GetX();
    u = static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ()));
    v = static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ()));
    w = static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ()));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &DlVertexingAlgorithm::GetTrueVertex() const
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

#ifdef MONITORING
void DlVertexingAlgorithm::PopulateRootTree(const std::vector<VertexTuple> &vertexTuples, const pandora::CartesianPointVector &vertexCandidatesU,
    const pandora::CartesianPointVector &vertexCandidatesV, const pandora::CartesianPointVector &vertexCandidatesW) const
{
    if (m_writeTree)
    {
        const MCParticleList *pMCParticleList{nullptr};
        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetCurrentList(*this, pMCParticleList))
        {
            if (pMCParticleList)
            {
                LArMCParticleHelper::MCContributionMap mcToHitsMap;
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
                        {
                            const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
                            const CartesianVector &trueVertex{primaries.front()->GetVertex()};
                            if (LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, "dune_fd_hd"))
                            {
                                const CartesianVector &recoVertex{vertexTuples.front().GetPosition()};
                                const float tx{trueVertex.GetX()};
                                const float tu{static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ()))};
                                const float tv{static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ()))};
                                const float tw{static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ()))};
                                const float rx_u{vertexCandidatesU.front().GetX()};
                                const float ru{vertexCandidatesU.front().GetZ()};
                                const float rx_v{vertexCandidatesV.front().GetX()};
                                const float rv{vertexCandidatesV.front().GetZ()};
                                const float rx_w{vertexCandidatesW.front().GetX()};
                                const float rw{vertexCandidatesW.front().GetZ()};
                                const float dr_u{std::sqrt((rx_u - tx) * (rx_u - tx) + (ru - tu) * (ru - tu))};
                                const float dr_v{std::sqrt((rx_v - tx) * (rx_v - tx) + (rv - tv) * (rv - tv))};
                                const float dr_w{std::sqrt((rx_w - tx) * (rx_w - tx) + (rw - tw) * (rw - tw))};
                                const CartesianVector &dv{recoVertex - trueVertex};
                                const float dr{dv.GetMagnitude()};
                                const float dx{dv.GetX()}, dy{dv.GetY()}, dz{dv.GetZ()};
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "event", m_event));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "pass", m_pass));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr_u", dr_u));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr_v", dr_v));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr_w", dr_w));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dr", dr));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dx", dx));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dy", dy));
                                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "dz", dz));
                                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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

DlVertexingAlgorithm::Canvas::Canvas(const int width, const int height, const int colOffset, const int rowOffset, const float xMin,
    const float xMax, const float zMin, const float zMax) :
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
    for (int r = 0; r < m_height; ++r)
        m_canvas[r] = new float[m_width]{};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::Canvas::~Canvas()
{
    for (int r = 0; r < m_height; ++r)
        delete[] m_canvas[r];
    delete[] m_canvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::VertexTuple::VertexTuple(
    const CartesianVector &vertexU, const CartesianVector &vertexV, const CartesianVector &vertexW) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    (void)vertexW;
    m_pos = CartesianVector(vertexU.GetX(), vertexU.GetY(), vertexV.GetZ());
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::VertexTuple::VertexTuple(
    const CartesianVector &vertex1, const CartesianVector &vertex2, const HitType view1, const HitType view2) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    if (view1 == TPC_VIEW_U && view2 == TPC_VIEW_V)
    {
        m_pos = CartesianVector(vertex1.GetX(), vertex1.GetY(), vertex2.GetZ());
    }
    else if (view1 == TPC_VIEW_U && view2 == TPC_VIEW_W)
    {
        m_pos = CartesianVector(vertex1.GetX(), vertex1.GetY(), vertex2.GetZ());
    }
    else
    {
        m_pos = CartesianVector(vertex1.GetX(), vertex2.GetY(), vertex1.GetZ());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &DlVertexingAlgorithm::VertexTuple::GetPosition() const
{
    return m_pos;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float DlVertexingAlgorithm::VertexTuple::GetChi2() const
{
    return m_chi2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::string DlVertexingAlgorithm::VertexTuple::ToString() const
{
    const float x{m_pos.GetX()}, y{m_pos.GetY()}, z{m_pos.GetZ()};
    return "3D pos: (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")   X2 = " + std::to_string(m_chi2);
}

} // namespace lar_dl_content
