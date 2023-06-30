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
    m_imageSize{256},
    m_driftStep{0.5f},
    m_findSecondaries{false},
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
    {
        if (m_pass == 1)
        {
            return this->PrepareTrainingSample();
        }
        else
        {
            const VertexList *pVertexList(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pVertexList));
            if (pVertexList->empty())
                return STATUS_CODE_NOT_FOUND;
            for (size_t i = 0; i < pVertexList->size(); ++i)
                this->PrepareTrainingSample(static_cast<int>(i));
        }
    }
    else
    {
        if (m_pass == 1)
        {
            return this->Infer();
        }
        else
        {
            const VertexList *pVertexList(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pVertexList));
            if (pVertexList->empty())
                return STATUS_CODE_NOT_FOUND;
            for (size_t i = 0; i < pVertexList->size(); ++i)
                this->Infer(static_cast<int>(i));
        }
    }

    return STATUS_CODE_SUCCESS;
}

StatusCode DlVertexingAlgorithm::PrepareTrainingSample(int vertexIndex)
{
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMCToHitsMap(mcToHitsMap));
    MCParticleList hierarchy;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CompleteMCHierarchy(mcToHitsMap, hierarchy));

    // Get boundaries for hits and make x dimension common
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitList3D, pCaloHitList));
    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    float xMin{0.f}, xMax{0.f}, yMin{0.f}, yMax{0.f}, zMin{0.f}, zMax{0.f};
    this->GetHitRegion(*pCaloHitList, xMin, xMax, yMin, yMax, zMin, zMax, vertexIndex);
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    if (pCaloHitList2D->empty())
        return STATUS_CODE_SUCCESS;

    CartesianPointVector vertices;
    for (const MCParticle *mc : hierarchy)
    {
        if (LArMCParticleHelper::IsNeutrino(mc))
        {
            vertices.emplace_back(mc->GetVertex());

            const MCParticleList &primaries{mc->GetDaughterList()};
            for (const MCParticle *pPrimary : primaries)
            {
                if (mcToHitsMap[pPrimary].empty())
                    continue;

                CartesianVector vertex(0, 0, 0);
                if (this->GetClearSecondaryVertex(pPrimary, *pCaloHitList2D, vertex))
                    vertices.emplace_back(vertex);
            }
        }
    }
    if (vertices.empty())
        return STATUS_CODE_SUCCESS;

    const std::string trainingFilename{m_trainingOutputFile + ".csv"};
    unsigned long nHits{0};

    LArMvaHelper::MvaFeatureVector featureVector;
    featureVector.emplace_back(static_cast<double>(m_event));

    int nContained{0};
    for (const CartesianVector &vertex : vertices)
    {
        const double xVtx{vertex.GetX()}, yVtx{vertex.GetY()}, zVtx{vertex.GetZ()};
        // Only train on events where the vertex resides within the image - with a small tolerance
        if (xVtx > (xMin - 1.f) && xVtx < (xMax + 1.f) && yVtx > (yMin - 1.f) && yVtx < (yMax + 1.f) &&
            zVtx > (zMin - 1.f) && zVtx < (zMax + 1.f))
        {
            featureVector.emplace_back(xVtx);
            featureVector.emplace_back(yVtx);
            featureVector.emplace_back(zVtx);
            ++nContained;
        }
    }
    if (!nContained)
        return STATUS_CODE_SUCCESS;
    featureVector.insert(featureVector.begin() + 1, static_cast<double>(nContained));

    // Retain the hit region
    featureVector.emplace_back(xMin);
    featureVector.emplace_back(xMax);
    featureVector.emplace_back(yMin);
    featureVector.emplace_back(yMax);
    featureVector.emplace_back(zMin);
    featureVector.emplace_back(zMax);

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()}, y{pCaloHit->GetPositionVector().GetY()},
            z{pCaloHit->GetPositionVector().GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};
        // If on a refinement pass, drop hits outside the region of interest
        if (m_pass > 1 && (x < xMin || x > xMax || y < yMin || y > yMax || z < zMin || z > zMax))
            continue;
        featureVector.emplace_back(static_cast<double>(x));
        featureVector.emplace_back(static_cast<double>(y));
        featureVector.emplace_back(static_cast<double>(z));
        featureVector.emplace_back(static_cast<double>(adc));
        ++nHits;
    }
    featureVector.insert(featureVector.begin() + 8 + 3 * nContained, static_cast<double>(nHits));
    // Only write out the feature vector if there were enough hits in the region of interest
    if (nHits > 10)
        LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::Infer(const int vertexIndex)
{
    // Get boundaries for hits and make x dimension common
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitList3D, pCaloHitList));
    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    float xMin{0.f}, xMax{0.f}, yMin{0.f}, yMax{0.f}, zMin{0.f}, zMax{0.f};
    this->GetHitRegion(*pCaloHitList, xMin, xMax, yMin, yMax, zMin, zMax, vertexIndex);

    Canvas canvasXY(m_imageSize, xMin, xMax, yMin, yMax);
    Canvas canvasXZ(m_imageSize, xMin, xMax, zMin, zMax);
    Canvas canvasYZ(m_imageSize, yMin, yMax, zMin, zMax);

    CartesianPointVector vertexCandidatesXY, vertexCandidatesXZ, vertexCandidatesYZ;
    for (const Projection projection : {XY, XZ, YZ})
    {
        std::cout << "Projection: " << projection << std::endl;
        LArDLHelper::TorchInput input;
        PixelVector pixelVector;
        switch (projection)
        {
            case XY:
                canvasXY.MakeNetworkInput(m_pass, *pCaloHitList, projection, input, pixelVector);
                break;
            case XZ:
                canvasXZ.MakeNetworkInput(m_pass, *pCaloHitList, projection, input, pixelVector);
                break;
            default:
                canvasYZ.MakeNetworkInput(m_pass, *pCaloHitList, projection, input, pixelVector);
                break;
        }

        // Run the input through the trained model
        LArDLHelper::TorchInputVector inputs;
        inputs.push_back(input);
        LArDLHelper::TorchOutput output;
        switch (projection)
        {
            case XY:
                LArDLHelper::Forward(m_modelXY, inputs, output);
                break;
            case XZ:
                LArDLHelper::Forward(m_modelXZ, inputs, output);
                break;
            default:
                LArDLHelper::Forward(m_modelYZ, inputs, output);
                break;
        }

        // we want the maximum value in the num_classes dimension (1) for every pixel
        auto classes{torch::argmax(output, 1)};
        // the argmax result is a 1 x height x width tensor where each element is a class id
        auto classesAccessor{classes.accessor<long, 3>()};
        const double scaleFactor{std::sqrt(2 * m_imageSize * m_imageSize)};
        std::map<int, bool> haveSeenMap;
        for (const auto &[row, col] : pixelVector)
        {
            const auto cls{classesAccessor[0][row][col]};
            if (cls > 0 && cls < m_nClasses)
            {
                const int inner{static_cast<int>(std::round(std::ceil(scaleFactor * m_thresholds[cls - 1])))};
                const int outer{static_cast<int>(std::round(std::ceil(scaleFactor * m_thresholds[cls])))};
                if (inner < 0.8f)
                {
                    std::cout << "row: " << row << " col:" << col << " cls: " << cls << " i: " << inner << " o: "  << outer <<
                        " = " << (1.f / (outer * outer - inner * inner)) << std::endl;
                    switch (projection)
                    {
                        case XY:
                            canvasXY.DrawRing(row, col, inner, outer, 1.f / (outer * outer - inner * inner));
                            break;
                        case XZ:
                            canvasXZ.DrawRing(row, col, inner, outer, 1.f / (outer * outer - inner * inner));
                            break;
                        default:
                            canvasYZ.DrawRing(row, col, inner, outer, 1.f / (outer * outer - inner * inner));
                            break;
                    }
                }
            }
        }

        switch (projection)
        {
            case XY:
                canvasXY.ExtractCandidateVertices(vertexCandidatesXY);
                break;
            case XZ:
                canvasXZ.ExtractCandidateVertices(vertexCandidatesXY);
                break;
            default:
                canvasYZ.ExtractCandidateVertices(vertexCandidatesXY);
                break;
        }
    }

    // Process candidates here - this will need some reworking
    /*CartesianPointVector vertexCandidates;
    for (const VertexTuple &tuple : vertexTuples)
        vertexCandidates.emplace_back(tuple.GetPosition());
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->MakeCandidateVertexList(vertexCandidates));*/

    return STATUS_CODE_SUCCESS;
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

void DlVertexingAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &yMin, float &yMax, float &zMin, float &zMax,
    const int vertexIndex) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    yMin = std::numeric_limits<float>::max();
    yMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
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

    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    // Choice of wire pitch is arbitrary here, so choosing collection plane
    const float pitch(pTPC->GetWirePitchW());

    if (m_pass > 1)
    {
        const VertexList *pVertexList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pVertexList));
        if (pVertexList->empty())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        auto iter{pVertexList->begin()};
        std::advance(iter, vertexIndex);
        const CartesianVector &vertex{(*iter)->GetPosition()};

        // Get hit distribution left/right asymmetry
        this->GetAsymmetryBounds(caloHitList, vertex, xMin, xMax, yMin, yMax, zMin, zMax);
    }

    // Avoid unreasonable rescaling of very small hit regions, pixels are assumed to be 0.5cm in x and wire pitch in z
    // ATTN: Rescaling is to a size 1 pixel smaller than the intended image to ensure all hits fit within an imaged binned
    // to be one pixel wider than this
    const float xRange{xMax - xMin}, zRange{zMax - zMin};
    const float minXSpan{m_driftStep * (m_imageSize - 1)};
    if (xRange < minXSpan)
    {
        const float padding{0.5f * (minXSpan - xRange)};
        xMin -= padding;
        xMax += padding;
    }
    const float minZSpan{pitch * (m_imageSize - 1)};
    if (zRange < minZSpan)
    {
        const float padding{0.5f * (minZSpan - zRange)};
        zMin -= padding;
        zMax += padding;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::GetAsymmetryBounds(const CaloHitList &caloHitList, const CartesianVector &vertex, float &xMin, float &xMax,
    float &yMin, float &yMax, float &zMin, float &zMax) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    // Choice of wire pitch is arbitrary here, so choosing collection plane
    const float pitch(pTPC->GetWirePitchW());

    // Get hit distribution left/right asymmetry
    int nHitsLeft{0}, nHitsRight{0};
    const double xVtx{vertex.GetX()};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        if (pos.GetX() <= xVtx)
            ++nHitsLeft;
        else
            ++nHitsRight;
    }
    int nHitsTotal{nHitsLeft + nHitsRight};
    if (nHitsTotal == 0)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    const float xAsymmetry{nHitsLeft / static_cast<float>(nHitsTotal)};

    // Get hit distribution up/down asymmetry
    int nHitsUp{0}, nHitsDown{0};
    const double yVtx{vertex.GetY()};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        if (pos.GetY() <= yVtx)
            ++nHitsDown;
        else
            ++nHitsUp;
    }
    nHitsTotal = nHitsUp + nHitsDown;
    if (nHitsTotal == 0)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    const float yAsymmetry{nHitsDown / static_cast<float>(nHitsTotal)};

    // Get hit distribution upstream/downstream asymmetry
    int nHitsUpstream{0}, nHitsDownstream{0};
    const double zVtx{vertex.GetZ()};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        if (pos.GetZ() <= zVtx)
            ++nHitsUpstream;
        else
            ++nHitsDownstream;
    }
    nHitsTotal = nHitsUpstream + nHitsDownstream;
    if (nHitsTotal == 0)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    const float zAsymmetry{nHitsUpstream / static_cast<float>(nHitsTotal)};

    // width/height choice here is a little arbitrary because of the 3 orthogonal views, but essentially assumes a square image
    const float xSpan{m_driftStep * (m_imageSize - 1)};
    xMin = xVtx - xAsymmetry * xSpan;
    xMax = xMin + xSpan;
    const float ySpan{pitch * (m_imageSize - 1)};
    yMin = yVtx - yAsymmetry * ySpan;
    yMax = yMin + ySpan;
    const float zSpan{pitch * (m_imageSize - 1)};
    zMin = zVtx - zAsymmetry * zSpan;
    zMax = zMin + zSpan;
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

bool DlVertexingAlgorithm::GetClearSecondaryVertex(const MCParticle *const pParent, const CaloHitList &caloHitList, CartesianVector &vertex) const
{
    const MCParticleList &secondaries{pParent->GetDaughterList()};
    const MCParticle *pLastSecondary{nullptr};
    int numVisible{0};
    for (const MCParticle *pSecondary : secondaries)
    {
        for (const CaloHit *pCaloHit : caloHitList)
        {
            try
            {
                const MCParticle *pMatchedMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (pSecondary == pMatchedMC)
                {
                    pLastSecondary = pSecondary;
                    ++numVisible;
                    break;
                }
            }
            catch (...)
            {
            }
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
                return this->GetClearSecondaryVertex(pLastSecondary, caloHitList, vertex);
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

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Pass", m_pass));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageSize", m_imageSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FindSecondaries", m_findSecondaries));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DistanceThresholds", m_thresholds));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitList3D", m_caloHitList3D));
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
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameXY", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelXY);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameXZ", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelXZ);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameYZ", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelYZ);
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

DlVertexingAlgorithm::Canvas::Canvas(const int imageSize, const float minCol, const float maxCol, const float minRow, const float maxRow) :
    m_imageSize{imageSize},
    m_minCol{minCol},
    m_maxCol{maxCol},
    m_minRow{minRow},
    m_maxRow{maxRow},
    m_canvas{nullptr}
{
    m_canvas = new float *[m_imageSize];
    for (int row = 0; row < m_imageSize; ++row)
        m_canvas[row] = new float[m_imageSize]{};
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::Canvas::~Canvas()
{
    for (int row = 0; row < m_imageSize; ++row)
        delete[] m_canvas[row];
    delete[] m_canvas;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::Canvas::MakeNetworkInput(const int pass, const CaloHitList &caloHits, const Projection projection,
    LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const
{
    const float pitch(0.4790000021565f);
    const float driftStep{0.5f};

    // Determine the bin edges
    std::vector<double> xBinEdges(m_imageSize + 1);
    std::vector<double> zBinEdges(m_imageSize + 1);
    xBinEdges[0] = m_minCol - 0.5f * driftStep;
    const double dx = ((m_maxCol + 0.5f * driftStep) - xBinEdges[0]) / m_imageSize;
    std::cout << "dx: " << dx << std::endl;
    for (int i = 1; i < m_imageSize + 1; ++i)
    {
        xBinEdges[i] = xBinEdges[i - 1] + dx;
        std::cout << "Bin Edge " << i << ": " << xBinEdges[i] << std::endl;
    }
    zBinEdges[0] = m_minRow - 0.5f * pitch;
    const double dz = ((m_maxRow + 0.5f * pitch) - zBinEdges[0]) / m_imageSize;
    std::cout << "dz: " << dx << std::endl;
    for (int i = 1; i < m_imageSize + 1; ++i)
    {
        zBinEdges[i] = zBinEdges[i - 1] + dz;
        std::cout << "Bin Edge " << i << ": " << zBinEdges[i] << std::endl;
    }

    LArDLHelper::InitialiseInput({1, 1, m_imageSize, m_imageSize}, networkInput);
    auto accessor = networkInput.accessor<float, 4>();

    for (const CaloHit *pCaloHit : caloHits)
    {
        const float x{projection == YZ ? pCaloHit->GetPositionVector().GetY() : pCaloHit->GetPositionVector().GetX()};
        const float z{projection == XY ? pCaloHit->GetPositionVector().GetY() : pCaloHit->GetPositionVector().GetZ()};
        std::cout << "Calo " << pCaloHit->GetPositionVector().GetX() << " " << pCaloHit->GetPositionVector().GetY() << " " <<
            pCaloHit->GetPositionVector().GetZ() << std::endl;
        std::cout << "   Proj: " << x << " " << z << " adc " << pCaloHit->GetMipEquivalentEnergy() << std::endl;
        if (pass > 1)
        {
            if (x < m_minCol || x > m_maxCol || z < m_minRow || z > m_maxRow)
                continue;
        }
        const float adc{pCaloHit->GetMipEquivalentEnergy()};
        const int pixelX{static_cast<int>(std::floor((x - xBinEdges[0]) / dx))};
        const int pixelZ{static_cast<int>(std::floor((z - zBinEdges[0]) / dz))};
        std::cout << "Pixel: " << pixelX << " " << pixelZ << std::endl;
        accessor[0][0][pixelZ][pixelX] += adc;
    }
    for (int row = 0; row < m_imageSize; ++row)
    {
        for (int col = 0; col < m_imageSize; ++col)
        {
            const float value{accessor[0][0][row][col]};
            if (value > 0)
                pixelVector.emplace_back(std::make_pair(row, col));
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::Canvas::DrawRing(const int row, const int col, const int inner, const int outer, const float weight) const
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
            int a1{row + rp2}, b1{col + c}, a2{row - rp2}, b2{col - c};
            int a3{row + c}, b3{col + rp2}, a4{row - c}, b4{col - rp2};
            if (0 <= a1 && a1 < m_imageSize && 0 <= b1 && b1 < m_imageSize)
                m_canvas[a1][b1] += weight;
            if (rp2 != c && 0 <= a3 && a3 < m_imageSize && 0 <= b3 && b3 < m_imageSize)
                m_canvas[a3][b3] += weight;
            if (rp2 != 0 && cp2 != 0)
            {
                if (0 <= a2 && a2 < m_imageSize && 0 <= b2 && b2 < m_imageSize)
                    m_canvas[a2][b2] += weight;
                if (rp2 != c && 0 <= a4 && a4 < m_imageSize && 0 <= b4 && b4 < m_imageSize)
                    m_canvas[a4][b4] += weight;
            }
            if (rp2 != 0)
            {
                if (0 <= a2 && a2 < m_imageSize && 0 <= b1 && b1 < m_imageSize)
                    m_canvas[a2][b1] += weight;
                if (rp2 != c && 0 <= a3 && a3 < m_imageSize && 0 <= b4 && b4 < m_imageSize)
                    m_canvas[a3][b4] += weight;
            }
            if (cp2 != 0)
            {
                if (0 <= a1 && a1 < m_imageSize && 0 <= b2 && b2 < m_imageSize)
                    m_canvas[a1][b2] += weight;
                if (rp2 != c && 0 <= a4 && a4 < m_imageSize && 0 <= b3 && b3 < m_imageSize)
                    m_canvas[a4][b3] += weight;
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

void DlVertexingAlgorithm::Canvas::ExtractCandidateVertices(pandora::CartesianPointVector &candidateVertices) const
{
    std::set<float> values;
    for (int row = 0; row < m_imageSize; ++row)
    {
        for (int col = 0; col < m_imageSize; ++col)
        {
            if (m_canvas[row][col] > 0)
            {
                std::cout << "(" << row << "," << col << "): " << m_canvas[row][col] << std::endl;
                values.insert(m_canvas[row][col]);
            }
        }
    }

    std::cout << "Values:" << std::endl;
    for (const float value : values)
        std::cout << value << std::endl;
    std::cout << std::endl;

    (void)candidateVertices;

/*    for (int row = 0; row < canvas.m_imageHeight; ++row)
    {
        for (int col = 0; col < canvas.m_imageWidth; ++col)
        {
            if (heatMap[row][col] >= 0.25f * best)
            {
                bool localMaximum{true};
                for (int r = -1; r <= 1; ++r)
                {
                    if (0 <= (row + r) && (row + r) < canvas.m_imageHeight)
                    {
                        for (int c = -1; c <= 1; ++c)
                        {
                            if ((r == 0) && (c == 0))
                                continue;
                            if (0 <= (col + c) && (col + c) < canvas.m_imageWidth)
                            {
                                if (heatMap[row][col] < heatMap[row + r][col + c])
                                {
                                    localMaximum = false;
                                    break;
                                }
                            }
                        }
                        if (!localMaximum)
                            break;
                    }
                }
                if (localMaximum)
                {
                    const float x0{static_cast<float>(col * dx + canvas.m_xMin)};
                    const float z0{static_cast<float>(row * dz + canvas.m_zMin)};
                    CartesianVector pt0(x0, 0.f, z0);

                    bool canAdd{true};
                    CartesianPointVector removals;
                    for (const CartesianVector &existing : positionVector)
                    {
                        if (pt0.GetDistanceSquared(existing) < 25.f)
                        {
                            const int row1{static_cast<int>((existing.GetZ() - canvas.m_zMin) / dz)};
                            const int col1{static_cast<int>((existing.GetX() - canvas.m_xMin) / dx)};

                            // Check heat value
                            if (heatMap[row][col] <= heatMap[row1][col1])
                            {
                                canAdd = false;
                                break;
                            }
                            else
                            {
                                removals.emplace_back(existing);
                            }
                        }
                    }
                    if (canAdd)
                    {
                        for (const CartesianVector &r : removals)
                        {
                            auto iter{std::find(positionVector.begin(), positionVector.end(), r)};
                            if (iter != positionVector.end())
                                positionVector.erase(iter);
                        }
                        positionVector.emplace_back(pt0);
                    }
                }
            }
        }
    }*/
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexingAlgorithm::Canvas::Update(const int radius2, int &col, int &row) const
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
//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::VertexTuple::VertexTuple(
    const Pandora &pandora, const CartesianVector &vertexU, const CartesianVector &vertexV, const CartesianVector &vertexW) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f},
    m_nInputs{3}
{
    m_inputs[TPC_VIEW_U] = &vertexU;
    m_inputs[TPC_VIEW_V] = &vertexV;
    m_inputs[TPC_VIEW_W] = &vertexW;

    LArGeometryHelper::MergeThreePositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vertexU, vertexV, vertexW, m_pos, m_chi2);
    if (m_chi2 > 1.f)
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
        }
        if (chi2UW < m_chi2)
        {
            m_pos = vertexUW;
            m_chi2 = chi2UW;
        }
        if (chi2VW < m_chi2)
        {
            m_pos = vertexVW;
            m_chi2 = chi2VW;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlVertexingAlgorithm::VertexTuple::VertexTuple(
    const Pandora &pandora, const CartesianVector &vertex1, const CartesianVector &vertex2, const HitType view1, const HitType view2) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f},
    m_nInputs{2}
{
    m_inputs[view1] = &vertex1;
    m_inputs[view2] = &vertex2;
    LArGeometryHelper::MergeTwoPositions3D(pandora, view1, view2, vertex1, vertex2, m_pos, m_chi2);
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

bool DlVertexingAlgorithm::VertexTuple::IsInput(const pandora::HitType view, const pandora::CartesianVector &vertex) const
{
    auto iter{m_inputs.find(view)};
    if (iter != m_inputs.end())
        return *(m_inputs.at(view)) == vertex;
    else
        return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

int DlVertexingAlgorithm::VertexTuple::GetNumInputs() const
{
    return m_nInputs;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float DlVertexingAlgorithm::VertexTuple::GetSeparation(const VertexTuple &other) const
{
    return (m_pos - other.GetPosition()).GetMagnitude();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::string DlVertexingAlgorithm::VertexTuple::ToString() const
{
    const float x{m_pos.GetX()}, y{m_pos.GetY()}, z{m_pos.GetZ()};
    return "3D pos [" + std::to_string(m_nInputs) + "]: (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) +
        ")   X2 = " + std::to_string(m_chi2);
}

} // namespace lar_dl_content
