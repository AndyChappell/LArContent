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
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
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
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    MCParticleList primariesList(primaries.begin(), primaries.end());

    InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primariesList)};
    if (!descriptor.IsDIS())
        throw StatusCodeException(STATUS_CODE_FAILURE);
    std::cout << "Is DIS" << std::endl;

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
    int hitBins[nBins]{};
    float closestApproach[nBins]{};
    CaloHitVector caloHits[nBins];
    for (int bin = 0; bin < nBins; ++bin)
        closestApproach[bin] = std::numeric_limits<float>::max();
    const float stepAngle{2 * M_PI * binAngle / 360.f};
    const float phase{0};
    for (const CaloHit *const pCaloHit : selectedHits)
    {
        const CartesianVector delta{pCaloHit->GetPositionVector() - pos};
        float dr{delta.GetMagnitudeSquared()};
        float theta{std::atan2(delta.GetZ(), delta.GetX())};
        if (theta < 0)
            theta += 2 * static_cast<float>(M_PI);
        theta -= phase;
        if (theta < 0)
            theta = 2 * static_cast<float>(M_PI) - theta;

        const int bin{static_cast<int>(std::floor(theta / stepAngle))};
        ++hitBins[bin];
        if (dr < closestApproach[bin])
            closestApproach[bin] = dr;
        caloHits[bin].emplace_back(pCaloHit);
    }
    IntVector sortedBins;
    for (int bin = 0; bin < nBins; ++bin)
    {
        if (hitBins[bin] > 0)
            sortedBins.emplace_back(bin);
    }
    auto SortBinsFunc = [&hitBins] (const int bin1, const int bin2) { return hitBins[bin1] > hitBins[bin2]; };
    std::sort(sortedBins.begin(), sortedBins.end(), SortBinsFunc);
    auto SortCaloHitsFunc = [&pos] (const CaloHit *pCaloHit1, const CaloHit *pCaloHit2)
    {
        const float dr1{(pCaloHit1->GetPositionVector() - pos).GetMagnitudeSquared()};
        const float dr2{(pCaloHit2->GetPositionVector() - pos).GetMagnitudeSquared()};

        return dr1 < dr2;
    };
    for (CaloHitVector caloHitVector : caloHits)
    {
        std::sort(caloHitVector.begin(), caloHitVector.end(), SortCaloHitsFunc);
        std::cout << "Vector:" << std::endl;
        for (const CaloHit *const pCaloHit : caloHitVector)
        {
            const float dr{(pCaloHit->GetPositionVector() - pos).GetMagnitudeSquared()};
            std::cout << "   " << dr << std::endl;
        }
    }

    bool clusterMade{true};
    while (clusterMade)
    {
        clusterMade = false;
        for (const int bin : sortedBins)
        {
            int target{bin};
            int binAlt1{bin > 0 ? bin - 1 : nBins - 1};
            int binAlt2{bin < (nBins - 1) ? bin + 1 : 0};
            float delta{closestApproach[bin]};
            float deltaAlt{closestApproach[binAlt1]};
            if (deltaAlt < delta)
            {
                target = binAlt1;
                delta = deltaAlt;
            }
            deltaAlt = closestApproach[binAlt2];
            if (deltaAlt < delta)
            {
                target = binAlt2;
                delta = deltaAlt;
            }
            //clusterMade = call a function to do the work on the target hit vector
            if (clusterMade)
            {
                hitBins[target] = 0;
                std::sort(sortedBins.begin(), sortedBins.end(), SortBinsFunc);
            }
        }
    }
/*    for (int bin = 0; bin < nBins; ++bin)
    {
        const int count{hitBins[bin]};
        const float angle{bin * stepAngle + phase};
        const float approach{closestApproach[bin]};
 
        std::cout << bin << ": " << count << " " << approach << " " << angle << std::endl;
    }
    std::cout << "Done" << std::endl;*/
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &selectedHits, "hits", RED));

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