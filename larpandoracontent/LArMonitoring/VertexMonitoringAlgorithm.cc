/**
 *  @file   larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

VertexMonitoringAlgorithm::VertexMonitoringAlgorithm() :
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexMonitoringAlgorithm::~VertexMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(
        SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));

    this->VisualizeVertices();

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::VisualizeVertices() const
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));
   
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
    const ParticleFlowObject *pRecoNeutrino{nullptr};
    if (!primaries.empty())
    {
        const MCParticle *primary{primaries.front()};
        const MCParticleList &parents{primary->GetParentList()};
        if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
            pTrueNeutrino = parents.front();
    }

    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            pRecoNeutrino = pPfo;
            break;
        }
    }

    if (pRecoNeutrino && pTrueNeutrino)
    {
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
        const CartesianVector tu(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ())));
        const CartesianVector tv(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ())));
        const CartesianVector tw(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ())));

        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pRecoNeutrino)->GetPosition()};
        const CartesianVector ru(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(recoVertex.GetY(), recoVertex.GetZ())));
        const CartesianVector rv(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(recoVertex.GetY(), recoVertex.GetZ())));
        const CartesianVector rw(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(recoVertex.GetY(), recoVertex.GetZ())));

        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tu, "U true vertex", BLUE, 2));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tv, "V true vertex", BLUE, 2));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tw, "W true vertex", BLUE, 2));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &ru, "U reco vertex", RED, 2));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rv, "V reco vertex", RED, 2));
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rw, "W reco vertex", RED, 2));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

