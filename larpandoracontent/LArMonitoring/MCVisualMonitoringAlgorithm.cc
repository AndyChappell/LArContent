/**
 *  @file   larpandoracontent/LArMonitoring/MCVisualMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the MC visualisation monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/MCVisualMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

MCVisualMonitoringAlgorithm::MCVisualMonitoringAlgorithm() :
    m_foldToPrimaries(false),
    m_foldToShowers(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCVisualMonitoringAlgorithm::~MCVisualMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCVisualMonitoringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 1;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_minHitSharingFraction = 0;
    parameters.m_foldBackHierarchy = m_foldToPrimaries;
    parameters.m_foldToLeadingShower = m_foldToShowers;

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    std::map<std::string, CaloHitList> uHits, vHits, wHits;
    for (const auto [ key, value ] : keys)
    {
        uHits[value] = CaloHitList();
        vHits[value] = CaloHitList();
        wHits[value] = CaloHitList();
    }
    uHits["other"] = CaloHitList();
    vHits["other"] = CaloHitList();
    wHits["other"] = CaloHitList();

    for (const auto [ pMC, pCaloHits ] : targetMCParticleToHitsMap)
    {
        for (const CaloHit *pCaloHit : pCaloHits)
        {
            const HitType view{pCaloHit->GetHitType()};

            try
            {
                const int pdg{std::abs(pMC->GetParticleId())};
                std::string key("other");
                if (keys.find(pdg) != keys.end())
                    key = keys[pdg];

                if (view == HitType::TPC_VIEW_U)
                    uHits[key].push_back(pCaloHit);
                else if (view == HitType::TPC_VIEW_V)
                    vHits[key].push_back(pCaloHit);
                else
                    wHits[key].push_back(pCaloHit);
            }
            catch (const StatusCodeException&)
            {
                continue;
            }
        }
    }

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    std::map<std::string, Color> colors = {{"mu", PINK}, {"e", GREEN}, {"gamma", ORANGE},
        {"kaon", GRAY}, {"pi", CYAN}, {"p", BLUE}, {"other", BLACK}};

    for (const auto [ key, value ] : keys)
    {
        if (uHits[value].size() > 0)
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits[value], "u_" + value, colors[value]));
        if (vHits[value].size() > 0)
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits[value], "v_" + value, colors[value]));
        if (wHits[value].size() > 0)
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits[value], "w_" + value, colors[value]));
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCVisualMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToShowers", m_foldToShowers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

