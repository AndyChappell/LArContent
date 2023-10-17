/**
 *  @file   larpandoracontent/LArVertex/VertexTuningAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArVertex/VertexTuningAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

VertexTuningAlgorithm::VertexTuningAlgorithm() :
    m_caloHitListName{""},
    m_vertexListName{""},
    m_outputVertexListName{""}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexTuningAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const VertexList *pVertexList{nullptr};
    StatusCode status{PandoraContentApi::GetList(*this, m_vertexListName, pVertexList)};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, status);

    if (status == STATUS_CODE_NOT_INITIALIZED)
        return STATUS_CODE_SUCCESS;

    if (!pCaloHitList || pCaloHitList->empty() || !pVertexList || pVertexList->empty())
        return STATUS_CODE_SUCCESS;

    CartesianPointVector refinedVertices;
    for (const Vertex *const pVertex : *pVertexList)
    {
        CartesianVector vertex(pVertex->GetPosition());
        this->Refine(*pCaloHitList, vertex);
        refinedVertices.emplace_back(vertex);
    }

    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

    for (const CartesianVector &position : refinedVertices)
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

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexTuningAlgorithm::Refine(const CaloHitList &caloHitList, CartesianVector &vertex) const
{
    if (caloHitList.empty())
        return;

    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    CartesianVector vertexU(vertex.GetX(), 0, static_cast<float>(transform->YZtoU(vertex.GetY(), vertex.GetZ())));
    CartesianVector vertexV(vertex.GetX(), 0, static_cast<float>(transform->YZtoV(vertex.GetY(), vertex.GetZ())));
    CartesianVector vertexW(vertex.GetX(), 0, static_cast<float>(transform->YZtoW(vertex.GetY(), vertex.GetZ())));

    // Gather the hits within a given radii of the vertex - tracks should be straighter nearby their interaction vertex
    CaloHitVector selectedHitsU, selectedHitsV, selectedHitsW;
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const HitType view{pCaloHit->GetHitType()};
        const CartesianVector &pos{view == TPC_VIEW_U ? vertexU : view == TPC_VIEW_V ? vertexV : vertexW};
        const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
        const float distanceSquared{hitPos.GetDistanceSquared(pos)};
        if (distanceSquared <= 9.f)
        {
            switch (view)
            {
                case TPC_VIEW_U:
                    selectedHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    selectedHitsV.emplace_back(pCaloHit);
                    break;
                default:
                    selectedHitsW.emplace_back(pCaloHit);
                    break;
            }
        }
    }
    std::sort(selectedHitsU.begin(), selectedHitsU.end(), LArClusterHelper::SortHitsByPositionInX);
    std::sort(selectedHitsV.begin(), selectedHitsV.end(), LArClusterHelper::SortHitsByPositionInX);
    std::sort(selectedHitsW.begin(), selectedHitsW.end(), LArClusterHelper::SortHitsByPositionInX);

    std::vector<std::tuple<const CartesianVector, const CartesianVector, const CartesianVector>> tuples;
    const CaloHit *pMatchV{nullptr}, *pMatchW{nullptr};
    CartesianVector posU(0, 0, 0), posV(0, 0, 0), posW(0, 0, 0);
    for (const CaloHit *const pCaloHitU : selectedHitsU)
    {
        float chi2Min{std::numeric_limits<float>::max()};
        const float xU{pCaloHitU->GetPositionVector().GetX()};
        const float xMinU{xU - 2 * pCaloHitU->GetCellSize1()};
        const float xMaxU{xU + 2 * pCaloHitU->GetCellSize1()};
        for (const CaloHit *const pCaloHitV : selectedHitsV)
        {
            const float xV{pCaloHitV->GetPositionVector().GetX()};
            if (xV < xMinU)
                continue;
            if (xV > xMaxU)
                break;
            for (const CaloHit *const pCaloHitW : selectedHitsW)
            {
                const float xW{pCaloHitW->GetPositionVector().GetX()};
                if (xW < xMinU)
                    continue;
                if (xW > xMaxU)
                    break;
                float chi2{0.f};
                LArGeometryHelper::MergeThreePositions(this->GetPandora(), pCaloHitU->GetPositionVector(), pCaloHitV->GetPositionVector(),
                    pCaloHitW->GetPositionVector(), posU, posV, posW, chi2);
                if (chi2 < chi2Min)
                {
                    chi2Min = chi2;
                    pMatchV = pCaloHitV;
                    pMatchW = pCaloHitW;
                }
            }
        }
        if (chi2Min < 1e-2)
        {
            LArGeometryHelper::MergeThreePositions(this->GetPandora(), pCaloHitU->GetPositionVector(), pMatchV->GetPositionVector(),
                pMatchW->GetPositionVector(), posU, posV, posW, chi2Min);
            tuples.emplace_back(std::make_tuple(posU, posV, posW));
        }
    }

    float maxScore{0.f};
    CartesianVector bestU(0, 0, 0), bestV(0, 0, 0), bestW(0, 0, 0);
    for (auto tuple : tuples)
    {
        selectedHitsU.clear();
        selectedHitsV.clear();
        selectedHitsW.clear();
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const HitType view{pCaloHit->GetHitType()};
            const CartesianVector &pos{view == TPC_VIEW_U ? std::get<0>(tuple) : view == TPC_VIEW_V ? std::get<1>(tuple) : std::get<2>(tuple)};
            const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
            const float distanceSquared{hitPos.GetDistanceSquared(pos)};
            if (distanceSquared <= 9.f)
            {
                switch (view)
                {
                    case TPC_VIEW_U:
                        selectedHitsU.emplace_back(pCaloHit);
                        break;
                    case TPC_VIEW_V:
                        selectedHitsV.emplace_back(pCaloHit);
                        break;
                    default:
                        selectedHitsW.emplace_back(pCaloHit);
                        break;
                }
            }
        }

        float score{this->GetScore(selectedHitsU, std::get<0>(tuple))};
        score += this->GetScore(selectedHitsV, std::get<1>(tuple));
        score += this->GetScore(selectedHitsW, std::get<2>(tuple));

        if (score > maxScore)
        {
            maxScore = score;
            bestU = std::get<0>(tuple);
            bestV = std::get<1>(tuple);
            bestW = std::get<2>(tuple);
        }
    }

    if (maxScore > 0)
    {
        float chi2{0.f};
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, bestU, bestV, bestW, vertex, chi2);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexTuningAlgorithm::GetScore(const CaloHitVector &caloHitVector, const CartesianVector &vertex) const
{
    double scale{0.};
    for (const CaloHit *const pCaloHit : caloHitVector)
        scale += pCaloHit->GetCellSize1();
    const double increment{1. / scale};

    const int binAngle{8};
    const int nBins{360 / binAngle};
    double hitBins[nBins]{};
    const double stepAngle{2 * M_PI * binAngle / 360.f};
    for (const CaloHit *const pCaloHit : caloHitVector)
    {
        const CartesianVector delta{pCaloHit->GetPositionVector() - vertex};
        double theta1{std::atan2(delta.GetZ(), delta.GetX())};
        if (theta1 < 0)
        {
            double theta2{theta1 + M_PI};
            theta1 += 2 * M_PI;
            const int bin1{static_cast<int>(std::floor(theta1 / stepAngle))};
            const int bin2{static_cast<int>(std::floor(theta2 / stepAngle))};
            hitBins[bin1] += increment * pCaloHit->GetCellSize1();
            hitBins[bin2] -= increment * pCaloHit->GetCellSize1();
        }
        else
        {
            double theta2{theta1 - M_PI};
            const int bin1{static_cast<int>(std::floor(theta1 / stepAngle))};
            const int bin2{static_cast<int>(std::floor(theta2 / stepAngle))};
            hitBins[bin1] += increment * pCaloHit->GetCellSize1();
            hitBins[bin2] -= increment * pCaloHit->GetCellSize1();
        }
    }

    double score{0.};
    for (int i = 0; i < nBins; ++i)
        score += hitBins[i] * hitBins[i];

    return static_cast<float>(score);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexTuningAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
