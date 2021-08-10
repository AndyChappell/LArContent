/**
 *  @file   larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

HierarchyValidationAlgorithm::HierarchyValidationAlgorithm() :
    m_writeTree{false},
    m_foldToPrimaries{false},
    m_foldToLeadingShowers{false},
    m_validateEvent{false},
    m_validateMC{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyValidationAlgorithm::~HierarchyValidationAlgorithm()
{
    if (m_writeTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    LArHierarchyHelper::MCHierarchy mcHierarchy;
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, m_foldToPrimaries, m_foldToLeadingShowers, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(*pPfoList, m_foldToPrimaries, m_foldToLeadingShowers, recoHierarchy);
    LArHierarchyHelper::MatchInfo matchInfo;
    LArHierarchyHelper::MatchHierarchies(mcHierarchy, recoHierarchy, matchInfo);
    matchInfo.Print(mcHierarchy);

    this->WritePfos(matchInfo);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::WritePfos(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    static int event{0};
    LArHierarchyHelper::MCMatchesVector matches;
    matches.insert(matches.end(), matchInfo.GetGoodMatches().begin(), matchInfo.GetGoodMatches().end());
    matches.insert(matches.end(), matchInfo.GetAboveThresholdMatches().begin(), matchInfo.GetAboveThresholdMatches().end());
    //matches.insert(matches.end(), matchInfo.GetSubThresholdMatches().begin(), matchInfo.GetSubThresholdMatches().end());

    for (const auto &match : matches)
    {
        std::cout << "MC " << match.GetMC()->GetParticleId() << " " << match.GetRecoMatches().size() << std::endl;
        for (const auto pNode : match.GetRecoMatches())
        {
            CaloHitList caloHitsU, caloHitsV, caloHitsW;
            FloatVector u_x, u_z, u_e, v_x, v_z, v_e, w_x, w_z, w_e;
            IntVector u_pi, v_pi, w_pi;
            const CaloHitList &nodeHits(pNode->GetCaloHits());
            int nPionHits{0};
            for (const CaloHit *pCaloHit : nodeHits)
            {
                const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                const bool isPrimaryPion{std::abs(pMCParticle->GetParticleId()) == 211 && LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 1};
                if (isPrimaryPion)
                    ++nPionHits;
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                        caloHitsU.emplace_back(pCaloHit);
                        u_x.emplace_back(pCaloHit->GetPositionVector().GetX());
                        u_z.emplace_back(pCaloHit->GetPositionVector().GetZ());
                        u_e.emplace_back(pCaloHit->GetInputEnergy());
                        u_pi.emplace_back(isPrimaryPion);
                        break;
                    case TPC_VIEW_V:
                        caloHitsV.emplace_back(pCaloHit);
                        v_x.emplace_back(pCaloHit->GetPositionVector().GetX());
                        v_z.emplace_back(pCaloHit->GetPositionVector().GetZ());
                        v_e.emplace_back(pCaloHit->GetInputEnergy());
                        v_pi.emplace_back(isPrimaryPion);
                        break;
                    case TPC_VIEW_W:
                        caloHitsW.emplace_back(pCaloHit);
                        w_x.emplace_back(pCaloHit->GetPositionVector().GetX());
                        w_z.emplace_back(pCaloHit->GetPositionVector().GetZ());
                        w_e.emplace_back(pCaloHit->GetInputEnergy());
                        w_pi.emplace_back(isPrimaryPion);
                        break;
                    default:
                        break;
                }
            }
/*            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsU, "U", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsV, "V", GREEN));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsW, "W", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
            if (nPionHits < 10)
                continue;

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", event));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "u_x", &u_x));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "u_z", &u_z));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "u_energies", &u_e));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "u_pi", &u_pi));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "v_x", &v_x));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "v_z", &v_z));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "v_energies", &v_e));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "v_pi", &v_pi));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "w_x", &w_x));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "w_z", &w_z));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "w_energies", &w_e));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "w_pi", &w_pi));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }
    ++event;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateEvent", m_validateEvent));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateMC", m_validateMC));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
        if (!(m_validateEvent || m_validateMC))
        {
            std::cout << "Error: WriteTree requested but no tree names found" << std::endl;
            return STATUS_CODE_NOT_FOUND;
        }
        else if (m_validateEvent && m_validateMC)
        {
            std::cout << "Error: Both event-level and MC-level validation requested simulataneously" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
