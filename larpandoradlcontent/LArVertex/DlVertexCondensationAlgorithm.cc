/**
 *  @file   larpandoradlcontent/LArVertex/DlVertexCondensationAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning vertexing algorithm.
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoradlcontent/LArVertex/DlVertexCondensationAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlVertexCondensationAlgorithm::DlVertexCondensationAlgorithm() :
    m_visualize{true},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

DlVertexCondensationAlgorithm::~DlVertexCondensationAlgorithm()
{
    if (m_trainingMode)
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

StatusCode DlVertexCondensationAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlVertexCondensationAlgorithm::PrepareTrainingSample()
{
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        if (!(view == TPC_VIEW_U || view == TPC_VIEW_V || view == TPC_VIEW_W))
            return STATUS_CODE_NOT_ALLOWED;

        LArMCParticleHelper::MCContributionMap mcToHitsMap;
        LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap, true);
        // Remove photon MC particles with neutron parents to clean up noise
        this->FilterNeutronInducedParticles(mcToHitsMap);
        MCVertexMap mcToVertexMap;
        for (const auto &[pMC, _] : mcToHitsMap)
            this->GetProjectedTrueVertices(pMC, view, mcToVertexMap);
        if (mcToVertexMap.empty())
            continue;

        // Find the closest hit to each vertex, along with the distance between the hit and the vertex
        MCVertexMap mcToMatchedVertexMap;
        this->MatchHitToVertex(mcToHitsMap, mcToVertexMap, mcToMatchedVertexMap);

        if (m_visualize)
            this->VisualizeByMC(mcToHitsMap, mcToVertexMap, mcToMatchedVertexMap);

        // Now consolidate any duplicate vertices and associate the relevant hits to the consolidated vertex
        VertexHitsMap vertexHitsMap;
        this->ConsolidateVertices(mcToMatchedVertexMap, mcToHitsMap, vertexHitsMap);

        this->MakeTrainingFile(listname, vertexHitsMap, mcToMatchedVertexMap, *pCaloHitList, mcToHitsMap);

        if (m_visualize)
            this->VisualizeByVertex(vertexHitsMap);

        // Need to construct training file
        // When allocating hits to condensation points, any hit without a corresponding particle in the
        // mcToHitsMap should be labelled as background - note, we'll need to find those in the full hit list
        // because currently they will just be absent from the relevant maps
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::GetProjectedTrueVertices(const MCParticle *const pMC, const HitType view, MCVertexMap &mcVertexMap) const
{
    const LArTransformationPlugin *pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    const CartesianVector mcVertex{pMC->GetParticleId() == PHOTON ? pMC->GetEndpoint() : pMC->GetVertex()};
    CartesianVector viewVertex(0, 0, 0);
    switch (view)
    {
        case TPC_VIEW_U:
            viewVertex.SetValues(mcVertex.GetX(), 0, pTransform->YZtoU(mcVertex.GetY(), mcVertex.GetZ()));
            break;
        case TPC_VIEW_V:
            viewVertex.SetValues(mcVertex.GetX(), 0, pTransform->YZtoV(mcVertex.GetY(), mcVertex.GetZ()));
            break;
        case TPC_VIEW_W:
            viewVertex.SetValues(mcVertex.GetX(), 0, pTransform->YZtoW(mcVertex.GetY(), mcVertex.GetZ()));
            break;
        default:
            break;
    }
    // Emplace avoids issues with the lack of a default constructor for CartesianVector
    mcVertexMap.emplace(pMC, viewVertex);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::MatchHitToVertex(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCVertexMap &mcVertexMap,
    MCVertexMap &mcToMatchedVertexMap) const
{
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        const CartesianPointVector vertices(1, mcVertexMap.at(pMC));
        Eigen::MatrixXf hitMat(caloHitList.size(), 2), vertMat(1, 2);
        LArEigenHelper::Vectorize(caloHitList, hitMat);
        LArEigenHelper::Vectorize({vertices}, vertMat);

        if (hitMat.cols() != vertMat.cols())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        Eigen::VectorXf hitSq{hitMat.rowwise().squaredNorm()};
        Eigen::VectorXf vertSq{vertMat.rowwise().squaredNorm()};

        // Compute all pairwise squared distances:
        // D2(i,j) = |h_i|^2 + |v_j|^2 - 2 * h_i.dot(v_j)
        Eigen::MatrixXf dot{hitMat * vertMat.transpose()};
        Eigen::MatrixXf d2{(-2.0f * dot).colwise() + hitSq}; // add hit norms per column
        d2 = d2.rowwise() + vertSq.transpose();              // add vertex norms per row

        Eigen::VectorXf col{d2.col(0)};
        Eigen::Index idx{0};
        col.minCoeff(&idx);

        const int index{static_cast<int>(idx)};
        // Emplace avoids issues with the lack of a default constructor for CartesianVector
        mcToMatchedVertexMap.emplace(pMC, caloHitVector[index]->GetPositionVector());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::ConsolidateVertices(const MCVertexMap &mcToMatchedVertexMap, const LArMCParticleHelper::MCContributionMap &mcToHitsMap, VertexHitsMap &vertexHitsMap) const
{
    for (const auto &[pMC, vertex] : mcToMatchedVertexMap)
    {
        if (vertexHitsMap.find(vertex) == vertexHitsMap.end())
            vertexHitsMap[vertex] = CaloHitList();
        const CaloHitList &caloHitList{mcToHitsMap.at(pMC)};
        vertexHitsMap[vertex].insert(vertexHitsMap[vertex].end(), caloHitList.begin(), caloHitList.end());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::FilterNeutronInducedParticles(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    for (auto it = mcToHitsMap.begin(); it != mcToHitsMap.end();)
    {
        const MCParticle *const pMC{it->first};
        bool foundNeutronParent{false};
        if (std::abs(pMC->GetParticleId()) == PHOTON)
        {
            const MCParticle *pCurrentMC{pMC};
            while (!pCurrentMC->GetParentList().empty())
            {
                pCurrentMC = pCurrentMC->GetParentList().front();
                if (std::abs(pCurrentMC->GetParticleId()) == NEUTRON)
                {
                    it = mcToHitsMap.erase(it);
                    foundNeutronParent = true;
                    break;
                }
            }
        }
        if (!foundNeutronParent)
            ++it;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::MakeTrainingFile(const std::string &treeName, const VertexHitsMap &vertexHitsMap,
    const MCVertexMap &mcToMatchedVertexMap, const CaloHitList &fullCaloHitList, const LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    (void)treeName;
    HitVertexLabelMap hitToVertexLabelMap;
    HitVertexWeightMap hitToVertexWeightMap;
    CondensationPointMap hitToCondensationPointMap;

    const int BACKGROUND{-1}, CP{1};
    // Initialise maps
    for (const CaloHit *const pCaloHit : fullCaloHitList)
    {
        hitToVertexLabelMap[pCaloHit] = IntVector();
        hitToVertexWeightMap[pCaloHit] = FloatVector();
        hitToCondensationPointMap[pCaloHit] = BACKGROUND;
    }

    // Label all hit vertices as condensation points
    for (const auto &[vertex, caloHitList] : vertexHitsMap)
    {
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
            // The vertex position is identically equal to the corresponding hit position by construction
            if (hitPosition == vertex)
                hitToCondensationPointMap[pCaloHit] = CP;
        }
    }

    // Associate hits to condensation points
    for (const CaloHit *const pCaloHit : fullCaloHitList)
    {
        const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
        // Allow for fuzzy matching of hits to vertices
        // Get all associated MC particles for the hit
        const MCParticleWeightMap &mcParticleWeightMap{pCaloHit->GetMCParticleWeightMap()};
        for (const auto &[pMC, _] : mcParticleWeightMap)
        {
            // If an MC particle was filtered out, we want its hits to remain as background and not be associated to a vertex
            if (mcToHitsMap.find(pMC) == mcToHitsMap.end())
                continue;
            if (mcToMatchedVertexMap.find(pMC) != mcToMatchedVertexMap.end())
            {
                // Get the associated vertex for the current MC particle of the current hit
                const CartesianVector &matchedVertex{mcToMatchedVertexMap.at(pMC)};
                // Get all of the hits associated with the vertex
                const auto vertexIt{vertexHitsMap.find(matchedVertex)};
                if (vertexIt != vertexHitsMap.end())
                {
                    // If the current hit is associated to the vertex, label it with the vertex index (aligns with the vertexHitsMap order)
                    const int vertexIndex{static_cast<int>(std::distance(vertexHitsMap.begin(), vertexIt))};
                    const IntVector &currentLabels(hitToVertexLabelMap.at(pCaloHit));
                    // Only emplace if not already present (different MC particles can have the same vertex)
                    if (std::find(currentLabels.begin(), currentLabels.end(), vertexIndex) == currentLabels.end())
                    {
                        hitToVertexLabelMap[pCaloHit].emplace_back(vertexIndex);
                        // Store the distance between the hit and the vertex as weight (to be normalized later)
                        const float distance{(hitPosition - matchedVertex).GetMagnitude()};
                        hitToVertexWeightMap[pCaloHit].emplace_back(distance);
                    }
                }
            }
        }
    }

    // Print out some info
    int i{0};
    for (const auto &[vertex, _] : vertexHitsMap)
    {
        std::cout << "Vertex " << i << " at (" << vertex.GetX() << ", " << vertex.GetZ() << ")" << std::endl;
        ++i;
    }
    std::cout << std::endl;
    for (const auto &[pCaloHit, labels] : hitToVertexLabelMap)
    {
        std::cout << "Hit [CP: " << hitToCondensationPointMap.at(pCaloHit) << "] (" << pCaloHit->GetPositionVector().GetX() << ", " <<
            pCaloHit->GetPositionVector().GetZ() << ") - Labels: ";
        for (const int label : labels)
            std::cout << label << " ";
        std::cout << " Weights: ";
        const FloatVector &weights{hitToVertexWeightMap.at(pCaloHit)};
        for (const float weight : weights)
            std::cout << weight << " ";
        std::cout << std::endl;
    }
}


//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::VisualizeByMC(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCVertexMap &mcToVertexMap,
    const MCVertexMap &mcToMatchedVertexMap) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));

    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        std::cout << "MC Particle ID: " << pMC->GetParticleId() << " Hits: " << caloHitList.size() << std::endl;
        const MCParticle *pCurrentMC{pMC};
        while (!pCurrentMC->GetParentList().empty())
        {
            pCurrentMC = pCurrentMC->GetParentList().front();
            size_t nHits{mcToHitsMap.find(pCurrentMC) != mcToHitsMap.end() ? mcToHitsMap.at(pCurrentMC).size() : 0};
            std::cout << "  (" << pCurrentMC->GetParticleId() << ": " << nHits << ") -> ";
        }
        std::cout << "[]" << std::endl;

        const CartesianVector matchedVertex{mcToMatchedVertexMap.at(pMC)};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &matchedVertex, "MatchedHit", RED, 1));
        const CartesianVector mcVertex{mcToVertexMap.at(pMC)};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &mcVertex, "Vertex", BLUE, 3));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, std::to_string(std::abs(pMC->GetParticleId())), AUTOITER));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::VisualizeByVertex(const VertexHitsMap &vertexHitsMap) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    for (const auto &[vertex, caloHitList] : vertexHitsMap)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertex, "ConsolidatedVertex", GREEN, 1));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "Hits", AUTOITER));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexCondensationAlgorithm::Infer()
{
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexCondensationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    if (!m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
