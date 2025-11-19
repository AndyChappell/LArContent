/**
 *  @file   larpandoracontent/LArCheating/ValidationAlgorithm.cc
 *
 *  @brief  Implementation of the validation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMetrics/ValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ValidationAlgorithm::ValidationAlgorithm() :
    m_caloHitListName("CaloHitList2D"),
    m_mcParticleListName("Input"),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D", "NeutrinoParticles3D"}),
    m_fileName("Validation.root"),
    m_treeName("tree"),
    m_minPurity(0.8f),
    m_minCompleteness(0.65f),
    m_minRecoHits(30),
    m_minRecoHitsPerView(10),
    m_minRecoGoodViews(2),
    m_removeRecoNeutrons(true),
    m_selectRecoHits(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ValidationAlgorithm::~ValidationAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EventTree", m_fileName.c_str(), "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "PFPTree", m_fileName.c_str(), "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "TrackTree", m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ValidationAlgorithm::Run()
{
    std::cout << this->GetPandora().GetRun() << std::endl;
    std::cout << this->GetPandora().GetSubrun() << std::endl;
    std::cout << this->GetPandora().GetEvent() << std::endl;


    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    PfoList pfoList;

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
            continue;

        pfoList.insert(pfoList.begin(), pPfoList->begin(), pPfoList->end());
    }

    if (pfoList.empty())
        return STATUS_CODE_SUCCESS;

    // Do matching
    LArHierarchyHelper::FoldingParameters foldParameters;
    foldParameters.m_foldToLeadingShowers = true;
    const LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria recoCriteria(
        m_minRecoHits, m_minRecoHitsPerView, m_minRecoGoodViews, m_removeRecoNeutrons);
    LArHierarchyHelper::MCHierarchy mcHierarchy(recoCriteria);
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(pfoList, foldParameters, recoHierarchy);
    const LArHierarchyHelper::QualityCuts quality(m_minPurity, m_minCompleteness, m_selectRecoHits);
    LArHierarchyHelper::MatchInfo matchInfo(mcHierarchy, recoHierarchy, quality);
    LArHierarchyHelper::MatchHierarchies(matchInfo);

    // Get target MC and BestMatch
    MCParticleList nuParticles;
    mcHierarchy.GetRootMCParticles(nuParticles);
    MCParticleVector nuParticlesVec;
    nuParticlesVec.insert(nuParticlesVec.begin(), nuParticles.begin(), nuParticles.end());
    std::sort(nuParticlesVec.begin(), nuParticlesVec.end(), LArMCParticleHelper::SortByMomentum);

    // Input entry for each neutrino hierarchy
    for (unsigned int i = 0; i < nuParticlesVec.size(); ++i)
    {
        MCParticleVector targetMC; PfoVector bestRecoMatch;

        // ATTN: I am pretty certain this vector is sorted, check with Andy
        LArHierarchyHelper::MCMatchesVector mcMatchesVec(matchInfo.GetMatches(nuParticlesVec.at(i)));

        for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
        {
            targetMC.push_back(mcMatches.GetMC()->GetMCParticles().front());

            // Get matches
            const int nMatches(mcMatches.GetRecoMatches().size());

            if ((nMatches == 0) || (!mcMatches.IsQuality(quality)))
            {
                bestRecoMatch.push_back(nullptr);
            }
            else
            {
                bestRecoMatch.push_back(mcMatches.GetRecoMatches().front()->GetRecoParticles().front());
            }
        }

        // Run tools
        m_pEventValidationTool->Run(this, nuParticlesVec.at(i), targetMC);
        m_pPFPValidationTool->Run(this, nuParticlesVec.at(i), targetMC, bestRecoMatch);
        m_pTrackValidationTool->Run(this, targetMC, bestRecoMatch);
        m_pShowerValidationTool->Run(this, targetMC, bestRecoMatch);
    }



    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCompleteness", m_minCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHits", m_minRecoHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHitsPerView", m_minRecoHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoGoodViews", m_minRecoGoodViews));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RemoveRecoNeutrons", m_removeRecoNeutrons));


    AlgorithmTool *pAlgorithmTool1(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "EventValidation", pAlgorithmTool1));
    m_pEventValidationTool = dynamic_cast<EventValidationTool *>(pAlgorithmTool1);

    AlgorithmTool *pAlgorithmTool2(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "TrackValidation", pAlgorithmTool2));
    m_pTrackValidationTool = dynamic_cast<TrackValidationTool *>(pAlgorithmTool2);

    AlgorithmTool *pAlgorithmTool3(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "PFPValidation", pAlgorithmTool3));
    m_pPFPValidationTool = dynamic_cast<PFPValidationTool *>(pAlgorithmTool3);

    AlgorithmTool *pAlgorithmTool4(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerValidation", pAlgorithmTool4));
    m_pShowerValidationTool = dynamic_cast<ShowerValidationTool *>(pAlgorithmTool4);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
