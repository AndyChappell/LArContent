/**
 *  @file   larpandoracontent/LArMonitoring/MCVisualMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the MC visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/MCVisualMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

MCVisualMonitoringAlgorithm::MCVisualMonitoringAlgorithm() :
    m_foldToPrimaries(false),
    m_visualisePdg(true),
    m_visualiseProcess(false)
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

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    this->MakeSelection(pMCParticleList, pCaloHitList, targetMCParticleToHitsMap);

    if (m_visualisePdg)
        this->VisualiseByPdgCode(targetMCParticleToHitsMap);

    if (m_visualiseProcess)
        this->VisualiseByProcess(targetMCParticleToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCVisualMonitoringAlgorithm::VisualiseByPdgCode(const LArMCParticleHelper::MCContributionMap &mcMap)
{
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

    for (const auto [ pMC, pCaloHits ] : mcMap)
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
    std::map<std::string, Color> colors = {{"mu", PINK}, {"e", GREEN}, {"gamma", ORANGE}, {"kaon", GRAY}, {"pi", CYAN}, {"p", BLUE},
        {"other", BLACK}};

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCVisualMonitoringAlgorithm::VisualiseByProcess(const LArMCParticleHelper::MCContributionMap &mcMap)
{
    std::map<MCProcess, const std::string> keys = {
        {MC_PROC_INCIDENT_NU, "incident_nu"}, {MC_PROC_UNKNOWN, "unknown"}, {MC_PROC_PRIMARY, "primary"}, {MC_PROC_COMPT, "compt"},
        {MC_PROC_PHOT, "phot"}, {MC_PROC_ANNIHIL, "annihil"}, {MC_PROC_E_IONI, "e_ioni"}, {MC_PROC_E_BREM, "e_brem"},
        {MC_PROC_CONV, "conv"}, {MC_PROC_MU_IONI, "mu_ioni"}, {MC_PROC_MU_MINUS_CAPTURE_AT_REST, "mu_capt_at_rest"},
        {MC_PROC_NEUTRON_INELASTIC, "n_inelastic"}, {MC_PROC_N_CAPTURE, "n_capture"}, {MC_PROC_HAD_ELASTIC, "had_elastic"},
        {MC_PROC_DECAY, "decay"}, {MC_PROC_COULOMB_SCAT, "coulomb_scat"}, {MC_PROC_UNUSED, "unused"}, {MC_PROC_MU_BREM, "mu_brem"},
        {MC_PROC_MU_PAIR_PROD, "mu_pair_prod"}, {MC_PROC_PHOTON_INELASTIC, "photon_inelastic"}, {MC_PROC_HAD_IONI, "had_ioni"},
        {MC_PROC_PROTON_INELASTIC, "p_inelastic"}, {MC_PROC_PI_PLUS_INELASTIC, "pi_plus_inelastic"},
        {MC_PROC_CHIPS_NUCLEAR_CAPTURE_AT_REST, "CHIPS_nuc_capt_at_rest"}, {MC_PROC_PI_MINUS_INELASTIC, "pi_minus_inelastic"}
    };

    std::map<std::string, CaloHitList> uHits, vHits, wHits;
    for (const auto [ key, value ] : keys)
    {
        uHits[value] = CaloHitList();
        vHits[value] = CaloHitList();
        wHits[value] = CaloHitList();
    }

    for (const auto [ pMC, pCaloHits ] : mcMap)
    {
        for (const CaloHit *pCaloHit : pCaloHits)
        {
            const HitType view{pCaloHit->GetHitType()};

            try
            {
                const LArMCParticle *pLArMC(dynamic_cast<const LArMCParticle *>(pMC));
                const MCProcess proc(pLArMC->GetProcess());
                std::string key(keys[proc]);

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
    std::map<std::string, Color> colors = {
        {"incident_nu", WHITE}, {"unknown", GRAY}, {"primary", BLACK}, {"compt", RED},
        {"phot", GREEN}, {"annihil", GREEN}, {"e_ioni", BLUE}, {"e_brem", MAGENTA},
        {"conv", BLUE}, {"mu_ioni", BLUE}, {"mu_capt_at_rest", CYAN},
        {"n_inelastic", RED}, {"n_capture", CYAN}, {"had_elastic", RED},
        {"decay", ORANGE}, {"coulomb_scat", RED}, {"unused", WHITE}, {"mu_brem", MAGENTA},
        {"mu_pair_prod", GREEN}, {"photon_inelastic", RED}, {"had_ioni", BLUE},
        {"p_inelastic", RED}, {"pi_plus_inelastic", RED},
        {"CHIPS_nuc_capt_at_rest", CYAN}, {"pi_minus_inelastic", RED}
    };

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCVisualMonitoringAlgorithm::MakeSelection(const MCParticleList *pMCList, const CaloHitList *pCaloHitList,
    LArMCParticleHelper::MCContributionMap &mcMap)
{
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 1;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_minHitSharingFraction = 0;
    parameters.m_foldBackHierarchy = m_foldToPrimaries;

    LArMCParticleHelper::SelectReconstructableMCParticles(pMCList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCVisualMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries",
        m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualisePDG",
        m_visualisePdg));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualiseProcess",
        m_visualiseProcess));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

