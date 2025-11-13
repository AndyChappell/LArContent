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
        CartesianPointVector vertices;
        this->GetProjectedTrueVertices(mcToHitsMap, view, vertices);
        if (vertices.empty())
            continue;

        // Find the closest hit to each vertex, along with the distance between the hit and the vertex
        Eigen::MatrixXf hitMat(pCaloHitList->size(), 2), vertMat(vertices.size(), 2);
        LArEigenHelper::Vectorize(*pCaloHitList, hitMat);
        LArEigenHelper::Vectorize(vertices, vertMat);
        IntVector closestHitIndices;
        FloatVector closestHitDistances;
        this->MatchHitToVertex(hitMat, vertMat, closestHitIndices, closestHitDistances);

        if (m_visualize)
        {
            CaloHitVector hitVector(pCaloHitList->begin(), pCaloHitList->end());
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));

            for (const auto &[pMC, caloHitList] : mcToHitsMap)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, std::to_string(std::abs(pMC->GetParticleId())), AUTOITER));
            }
            for (const CartesianVector &vertex : vertices)
            {
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &vertex, "Vertex", BLUE, 3));
            }
            for (int index : closestHitIndices)
            {
                const CartesianVector &hitPos{hitVector[index]->GetPositionVector()};
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &hitPos, "MatchedHit", RED, 1));
            }

            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        // Need to construct training file
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::GetProjectedTrueVertices(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const HitType view,
    CartesianPointVector &vertices) const
{
    const LArTransformationPlugin *pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
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
        vertices.emplace_back(viewVertex);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlVertexCondensationAlgorithm::MatchHitToVertex(const Eigen::MatrixXf &hitMat, const Eigen::MatrixXf &vertMat,
    IntVector &closestHitIndices, FloatVector &closestHitDistances) const
{
    const size_t nVerts = static_cast<size_t>(vertMat.rows());

    // Ensure shapes make sense
    if (hitMat.cols() != vertMat.cols()) throw std::runtime_error("hit/vert dims mismatch");

    // Use VectorXf (column vector) for squared norms
    Eigen::VectorXf hitSq  = hitMat.rowwise().squaredNorm();   // (N_h x 1)
    Eigen::VectorXf vertSq = vertMat.rowwise().squaredNorm();  // (N_v x 1)

    // Compute all pairwise squared distances:
    // D2(i,j) = |h_i|^2 + |v_j|^2 - 2 * h_i.dot(v_j)
    Eigen::MatrixXf dot = hitMat * vertMat.transpose();        // (N_h x N_v)
    Eigen::MatrixXf d2 = (-2.0f * dot).colwise() + hitSq;      // add hit norms per column
    d2 = d2.rowwise() + vertSq.transpose();                    // add vertex norms per row

    // Resize outputs
    closestHitIndices.resize(nVerts);
    closestHitDistances.resize(nVerts);

    // Visualize the MC particles and their true vertices here, then visualize the hit-based vertices

    // For each vertex (column j in d2), find index i with minimal d2(i,j)
    for (size_t j = 0; j < nVerts; ++j)
    {
        Eigen::VectorXf col = d2.col(static_cast<int>(j)); // N_h x 1
        Eigen::Index idx;
        float minD2 = col.minCoeff(&idx);

        closestHitIndices[j] = static_cast<int>(idx);
        closestHitDistances[j] = std::sqrt(std::max(minD2, 0.0f)); // guard tiny negatives

        std::cout << "Vertex " << j
                  << " (" << vertMat(int(j), 0) << ", " << vertMat(int(j), 1) << ") "
                  << "matched to hit " << idx
                  << " (" << hitMat(int(idx), 0) << ", " << hitMat(int(idx), 1) << ") "
                  << "distance " << closestHitDistances[j] << "\n";
    }
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
