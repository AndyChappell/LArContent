/**
 *  @file   larpandoracontent/LArMonitoring/VertexAssociatedHitMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VertexAssociatedHitMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

VertexAssociatedHitMonitoringAlgorithm::VertexAssociatedHitMonitoringAlgorithm() :
    m_caloHitListName{""},
    m_vertexListName{""},
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssociatedHitMonitoringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const VertexList *pVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));

    if (!pCaloHitList || pCaloHitList->empty() || !pVertexList || pVertexList->empty())
        return STATUS_CODE_SUCCESS;

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));

    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (const Vertex *const pVertex : *pVertexList)
    {
        this->IdentifyTrackStubs(*pCaloHitList, *pVertex);
    
        const CartesianVector &position{pVertex->GetPosition()};
        const CartesianVector w(position.GetX(), 0.f, static_cast<float>(transform->YZtoW(position.GetY(), position.GetZ())));

        CaloHitList selectedHits;
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
            const float distanceSquared{hitPosition.GetDistanceSquared(w)};
            if (distanceSquared <= 100.f)
                selectedHits.emplace_back(pCaloHit);
        }

        for (int r = 0; r < 360; r += 8)
        {
            CartesianVector start(w);
            CartesianVector end1(std::cos(2 * 3.14159 * r / 360.f), 0, std::sin(2 * 3.14159 * r / 360.f));
            CartesianVector end2(std::cos(2 * 3.14159 * (r + 8) / 360.f), 0, std::sin(2 * 3.14159 * (r + 8) / 360.f));
            end1 *= 10;
            end2 *= 10;
            end1 += start;
            end2 += start;
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &end1, "ref", BLACK, 1, 1));
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &end2, "ref", BLACK, 1, 1));
        }

        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &selectedHits, "hits", RED));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &w, "vertex", BLUE, 1));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexAssociatedHitMonitoringAlgorithm::IdentifyTrackStubs(const CaloHitList &caloHitList, const Vertex &vertex) const
{
    if (caloHitList.empty())
        return;

    HitType view{caloHitList.front()->GetHitType()};
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    const CartesianVector &pos3D{vertex.GetPosition()};
    CartesianVector pos(0, 0, 0);
    std::cout << "View " << view << std::endl;
    switch (view)
    {
        case TPC_VIEW_U:
            pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoU(pos3D.GetY(), pos3D.GetZ())));
            break;
        case TPC_VIEW_V:
            pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoV(pos3D.GetY(), pos3D.GetZ())));
            break;
        case TPC_VIEW_W:
            pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoW(pos3D.GetY(), pos3D.GetZ())));
            break;
        default:
            return;
    }

    // Gather the hits within a given radii of the vertex - tracks should be straighter nearby their interaction vertex
    CaloHitList selectedHits;
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
        const float distanceSquared{hitPos.GetDistanceSquared(pos)};
        if (distanceSquared <= 100.f)
            selectedHits.emplace_back(pCaloHit);
    }

    // Use a narow angular binning to try to isolate the track stubs
    // Do this twice with a phase offset to reduce sensitivity to tracks straddling bins
    const int binAngle{8};
    const int nBins{360 / binAngle};
    int hitBins1[nBins]{};
    int hitBins2[nBins]{};
    const float stepAngle{2 * M_PI * binAngle / 360.f};
    const float phase1{0}, phase2{stepAngle / 2};
    for (const CaloHit *const pCaloHit : selectedHits)
    {
        const CartesianVector delta{pCaloHit->GetPositionVector() - pos};
        float theta1{std::atan2(delta.GetZ(), delta.GetX()) + static_cast<float>(M_PI)};
        float theta2{theta1};
        theta1 -= phase1;
        theta2 -= phase2;
        if (theta1 < 0)
            theta1 = 2 * static_cast<float>(M_PI) - theta1;
        if (theta2 < 0)
            theta2 = 2 * static_cast<float>(M_PI) - theta2;

        const int bin1{static_cast<int>(std::floor(theta1 / stepAngle))};
        const int bin2{static_cast<int>(std::floor(theta2 / stepAngle))};
        ++hitBins1[bin1];
        ++hitBins2[bin2];
    }

    std::cout << "Num bins: " << nBins << std::endl;
    for (int bin = 0; bin < nBins; ++bin)
    {
        std::cout << bin << ": " << hitBins1[bin] << "   " << hitBins2[bin] << std::endl;
    }
    std::cout << "Done" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssociatedHitMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
