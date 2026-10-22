/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/SliceMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

SliceMonitoringAlgorithm::SliceMonitoringAlgorithm() :
    m_caloHitListName("CaloHitListW"),
    m_vertexListName("NeutrinoVertices3D")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SliceMonitoringAlgorithm::~SliceMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            mcToHitsMap[pMC].emplace_back(pCaloHit);
        }
        catch (StatusCodeException &statusCodeException)
        {
        }
    }

    LArMCParticleHelper::MCContributionMap nuToHitsMap;
    CaloHitList cosmicHitList;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const MCParticle *pParent{pMC};
        while (!pParent->GetParentList().empty())
        {
            pParent = pParent->GetParentList().front();
        }
        const int pdg{std::abs(pParent->GetParticleId())};
        if (pdg == 12 || pdg == 14 || pdg == 16)
        {
            nuToHitsMap[pParent].insert(nuToHitsMap[pParent].end(), caloHitList.begin(), caloHitList.end());
        }
        else
        {
            cosmicHitList.insert(cosmicHitList.end(), caloHitList.begin(), caloHitList.end());
        }
    }

    const MCParticle *pTrueNeutrino{nullptr};
    size_t maxHits{0};
    for (const auto &[pMC, caloHitList] : nuToHitsMap)
    {
        if (caloHitList.size() >= 10)
        {
            if (caloHitList.size() < maxHits)
                continue;
            maxHits = caloHitList.size();
            pTrueNeutrino = pMC;
        }
    }
    if (!pTrueNeutrino)
        return STATUS_CODE_SUCCESS;

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    //const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueVertex, "True vertex (SCE)", BLUE, 2));
    const VertexList *pVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    if (pVertexList && !pVertexList->empty())
    {
        const CartesianVector &recoVertex{pVertexList->front()->GetPosition()};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &recoVertex, "Reco vertex", RED, 2));
    }
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &nuToHitsMap[pTrueNeutrino], "True Neutrino Hits", BLACK));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &cosmicHitList, "Cosmic Hits", GRAY));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
