/**
 *  @file   larpandoracontent/LArVertex/VertexRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the vertex refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArVertex/VertexRefinementAlgorithm.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

VertexRefinementAlgorithm::VertexRefinementAlgorithm() :
    m_caloHitListName{"CaloHitList2D"},
    m_hitRadii(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexRefinementAlgorithm::Run()
{
    const VertexList *pInputVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    if (!pInputVertexList || pInputVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexRefinementAlgorithm: unable to find current vertex list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexRefinementAlgorithm: unable to find calo hit list \'" << m_caloHitListName << "\'" << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    const VertexList *pOutputVertexList{nullptr};
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputVertexList, temporaryListName));

    this->RefineVertices(*pInputVertexList, *pCaloHitList);

    if (!pOutputVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexRefinementAlgorithm::RefineVertices(const VertexList &vertexList, const CaloHitList &caloHitList) const
{
    CaloHitVector caloHitVectorU, caloHitVectorV, caloHitVectorW;
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                caloHitVectorU.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_V:
                caloHitVectorV.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_W:
                caloHitVectorW.emplace_back(pCaloHit);
                break;
            default:
                break;
        }
    }

    for (const Vertex *const pVertex : vertexList)
    {

        const CartesianVector originalPosition(pVertex->GetPosition());
        const CartesianVector &originalVtxU{LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_U)};
        //const CartesianVector &originalVtxV{LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_V)};
        //const CartesianVector &originalVtxW{LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_W)};

        CaloHitList nearbyHitListU;
        this->GetNearbyHits(caloHitVectorU, originalVtxU, nearbyHitListU);
        const CartesianVector vtxU(this->RefineVertexTwoD(nearbyHitListU));
        (void)vtxU;

/*        CaloHitList caloHitListU, caloHitListV, caloHitListW;
        // Collect calo hits around this vertex here

        const CartesianVector vtxU(this->RefineVertexTwoD(caloHitListU, originalVtxU));
        const CartesianVector vtxV(this->RefineVertexTwoD(caloHitListV, originalVtxV));
        const CartesianVector vtxW(this->RefineVertexTwoD(caloHitListW, originalVtxW));

        CartesianVector vtxUV(0.f, 0.f, 0.f), vtxUW(0.f, 0.f, 0.f), vtxVW(0.f, 0.f, 0.f), vtx3D(0.f, 0.f, 0.f), position3D(0.f, 0.f, 0.f);
        float chi2UV(0.f), chi2UW(0.f), chi2VW(0.f), chi23D(0.f), chi2(0.f);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, vtxU, vtxV, vtxUV, chi2UV);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, vtxU, vtxW, vtxUW, chi2UW);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vtxV, vtxW, vtxVW, chi2VW);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vtxU, vtxV, vtxW, vtx3D, chi23D);

        if (chi2UV < chi2UW && chi2UV < chi2VW && chi2UV < chi23D)
        {
            position3D = vtxUV;
            chi2 = chi2UV;
        }
        else if (chi2UW < chi2VW && chi2UW < chi23D)
        {
            position3D = vtxUW;
            chi2 = chi2UW;
        }
        else if (chi2VW < chi23D)
        {
            position3D = vtxVW;
            chi2 = chi2VW;
        }
        else
        {
            position3D = vtx3D;
            chi2 = chi23D;
        }

        if (chi2 > m_chiSquaredCut)
            position3D = originalPosition;

        if ((position3D - originalPosition).GetMagnitude() > m_distanceCut)
            position3D = originalPosition;

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position3D;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pNewVertex(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));*/
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<class T>
void VertexRefinementAlgorithm::Vectorize(const T &caloHitContainer, Eigen::MatrixXf &hitMatrix, Eigen::RowVectorXf &weightVector) const
{
    int i{0};
    for (const CaloHit *const pCaloHit : caloHitContainer)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        hitMatrix(i, 0) = pos.GetX();
        hitMatrix(i, 1) = pos.GetZ();
        weightVector(i) = std::pow(1.f / (pCaloHit->GetCellSize1() / 0.5f), 2);
        ++i;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexRefinementAlgorithm::GetNearbyHits(const CaloHitVector &hitVector, const CartesianVector &centroid, CaloHitList &nearbyHitList) const
{
    Eigen::MatrixXf hitMatrix(hitVector.size(), 2);
    Eigen::RowVectorXf dummy(hitVector.size());
    this->Vectorize(hitVector, hitMatrix, dummy);
    Eigen::RowVectorXf vertex(2);
    vertex << centroid.GetX(), centroid.GetZ();
    Eigen::MatrixXf norms((hitMatrix.rowwise() - vertex).array().pow(2).rowwise().sum());
    for (int r = 0; r < hitMatrix.rows(); ++r)
    {
        if (norms(r, 0) < m_hitRadii * m_hitRadii)
            nearbyHitList.emplace_back(hitVector.at(r));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector VertexRefinementAlgorithm::RefineVertexTwoD(const CaloHitList &caloHitList) const
{
    const int nBins{36};  // 72 //180
    const float pi{static_cast<float>(M_PI)};
    Eigen::RowVectorXf piVec{Eigen::RowVectorXf::Constant(caloHitList.size(), pi)};
    Eigen::RowVectorXf zeroVec{Eigen::RowVectorXf::Zero(caloHitList.size())};
    Eigen::RowVectorXf unitVec{Eigen::RowVectorXf::Constant(caloHitList.size(), 1.f)};
    Eigen::MatrixXf hitMatrix(caloHitList.size(), 2);
    Eigen::RowVectorXf weightVector(caloHitList.size());
    this->Vectorize(caloHitList, hitMatrix, weightVector);
    Eigen::RowVectorXi results{Eigen::RowVectorXi::Zero(hitMatrix.rows())};
    float best{0.f};
    for (int r = 0; r < hitMatrix.rows(); ++r)
    {
        Eigen::RowVectorXf row(2);
        row << hitMatrix(r, 0), hitMatrix(r, 1);
        // Compute dx, dz and angle between each hit and the candidate vertex hit
        Eigen::MatrixXf deltas(hitMatrix.rowwise() - row);
        Eigen::RowVectorXf norms(deltas.array().pow(2).rowwise().sum());
        Eigen::RowVectorXf invNorms(25.f / norms.array());
        Eigen::RowVectorXf phis(hitMatrix.rows());
        for (int i = 0; i < deltas.rows(); ++i)
            phis(i) = std::atan2(deltas(i, 1), deltas(i, 0));
        // Compute 180 degree rotation to suppress mid-track hits
        Eigen::RowVectorXf phisRot = (phis.array() >= 0).select(-piVec, piVec);
        Eigen::RowVectorXf rWeightVector = (norms.array() >= 25.f).select(invNorms, unitVec);
        phisRot += phis;

        phis = (phis.array() >= pi).select(zeroVec, phis);
        phisRot = (phisRot.array() >= pi).select(zeroVec, phisRot);
        Eigen::RowVectorXf counts{Eigen::RowVectorXf::Zero(nBins)};
        Eigen::RowVectorXf counts2{Eigen::RowVectorXf::Zero(nBins)};
        phis = phis.array() + pi;
        phis = phis.array() * nBins / (2 * pi);
        phisRot = phisRot.array() + pi;
        phisRot = phisRot.array() * nBins / (2 * pi);

        std::map<int, int> populated;
        for (int i = 0; i < deltas.rows(); ++i)
        {
            if (populated.find(static_cast<int>(phis(i))) != populated.end())
                ++populated[static_cast<int>(phis(i))];
            else
                populated[static_cast<int>(phis(i))] = 1;
            counts(static_cast<int>(phis(i))) += weightVector(i) * rWeightVector(i);
            counts2(static_cast<int>(phisRot(i))) += weightVector(i) * rWeightVector(i);
        }
        int populatedBins{0};
        for (const auto &[key, value] : populated)
            populatedBins += value;
        if (populatedBins == 0)
            populatedBins = 1;
        results(r) = counts.array().pow(2).sum();// / populatedBins;
        if (results(r) > best)
        {
            best = results(r);
            for (int i = 0; i < nBins; ++i)
            {
                std::cout << i << ": " << counts(i) << "   " << counts2(i) << std::endl;
            }
            std::cout << "Total: " << results(r) << std::endl;

            const CartesianVector centroid(hitMatrix(r, 0), 0, hitMatrix(r, 1));
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "near", BLACK));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &centroid, "vtx", RED, 3));
            for (int i = 0; i < nBins; ++i)
            {
                const float theta{2.f * i * pi / nBins};
                const CartesianVector b{CartesianVector(std::cos(theta), 0, std::sin(theta)) * 2 * m_hitRadii + centroid};
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &centroid, &b, "bin " + std::to_string(i), GRAY, 1, 1));
            }
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
    Eigen::Index index;
    results.maxCoeff(&index);

    const CartesianVector centroid(hitMatrix(index, 0), 0, hitMatrix(index, 1));
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &centroid, "vtx", RED, 1));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "near", BLUE));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return centroid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitRadii", m_hitRadii));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
