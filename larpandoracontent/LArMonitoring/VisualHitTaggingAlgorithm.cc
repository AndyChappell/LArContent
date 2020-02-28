/**
 *  @file   larpandoracontent/LArDeepLearning/VisualHitTaggingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VisualHitTaggingAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

VisualHitTaggingAlgorithm::VisualHitTaggingAlgorithm() :
    m_caloHitListNames{},
    m_tagPid{false},
    m_tagPrimary{false},
    m_retainNonReconstructable{true}
{
}

//------------------------------------------------------------------------------

VisualHitTaggingAlgorithm::~VisualHitTaggingAlgorithm()
{
}

//------------------------------------------------------------------------------

StatusCode VisualHitTaggingAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true,
                DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_minHitSharingFraction = 0.f;

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(
                *this, pMCParticleList));
    const CaloHitList *pCaloHitListAllViews(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(
                *this, "CaloHitList2D", pCaloHitListAllViews));

    LArMCParticleHelper::MCContributionMap primaryMCParticleToHitsMapAllViews;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitListAllViews,
        parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, primaryMCParticleToHitsMapAllViews);

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector(
            {primaryMCParticleToHitsMapAllViews}, mcPrimaryVector);
    std::cout << "Primary vector: " << mcPrimaryVector.size() << std::endl;

    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        std::string viewStr;
        if (isU) viewStr = "U/";
        else if (isV) viewStr = "V/";
        else viewStr = "W/";

        // Unfolded
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectUnfoldedReconstructableMCParticles(pMCParticleList, pCaloHitList,
                parameters, targetMCParticleToHitsMap);
        // Fold back to primary
        LArMCParticleHelper::MCContributionMap primaryMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList,
                parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, primaryMCParticleToHitsMap);

        CaloHitList nonReconstructablePrimaryHits;
        CaloHitList nonReconstructableHits;
        CaloHitList pHits;
        CaloHitList eHits;
        CaloHitList muHits;
        CaloHitList gammaHits;
        CaloHitList piHits;
        CaloHitList pi0Hits;
        CaloHitList kHits;
        CaloHitList nHits;
        CaloHitList otherHits;
        LArMCParticleHelper::MCContributionMap primaryHits;
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                if(m_tagPid)
                {
                    // Throw away non-reconstructable hits
                    if (targetMCParticleToHitsMap.find(pMCParticle) == targetMCParticleToHitsMap.end())
                    {
                        nonReconstructableHits.push_back(pCaloHit);
                        if (!m_retainNonReconstructable) continue;
                    }
                    int pdg = pMCParticle->GetParticleId();
                    switch (abs(pdg))
                    {
                        case 11:
                            eHits.push_back(pCaloHit);
                            break;
                        case 13:
                            muHits.push_back(pCaloHit);
                            break;
                        case 22:
                            gammaHits.push_back(pCaloHit);
                            break;
                        case 111:
                            pi0Hits.push_back(pCaloHit);
                            break;
                        case 211:
                            piHits.push_back(pCaloHit);
                            break;
                        case 311:
                        case 313:
                        case 321:
                        case 323:
                            kHits.push_back(pCaloHit);
                            break;
                        case 2212:
                            pHits.push_back(pCaloHit);
                            break;
                        case 2112:
                            nHits.push_back(pCaloHit);
                            break;
                        default:
                            std::cout << "Found PDG " << pdg << std::endl;
                            otherHits.push_back(pCaloHit);
                            break;
                    }

                }
                if(m_tagPrimary)
                {
                    if (primaryMCParticleToHitsMap.find(pMCParticle) == primaryMCParticleToHitsMap.end())
                    {
                        nonReconstructablePrimaryHits.push_back(pCaloHit);
                        if (!m_retainNonReconstructable) continue;
                    }
                    if (primaryHits.find(pMCParticle) == primaryHits.end())
                    {
                        primaryHits.insert(std::make_pair(pMCParticle, CaloHitList()));
                    }
                    primaryHits[pMCParticle].push_back(pCaloHit);
                }
            }
            catch (...)
            {
                continue;
            }
        }
        if (m_tagPid)
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &nonReconstructableHits, "NonReco", GRAY));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &pHits, viewStr + "Proton", BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &eHits, viewStr + "Electron", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &muHits, viewStr + "Muon", MAGENTA));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &gammaHits, viewStr + "Gamma", YELLOW));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &piHits, viewStr + "Pion", GREEN));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &pi0Hits, viewStr + "Pi0", LIGHTGREEN));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &kHits, viewStr + "Kaon", CYAN));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &nHits, viewStr + "Neutron", DARKBLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &otherHits, viewStr + "Other", BLACK));
        }
        if (m_tagPrimary)
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                        &nonReconstructablePrimaryHits, viewStr + "NonRecoPrimary", GRAY));
            int iPrimary = 0;
            Color colors[] = {BLACK,RED,GREEN,BLUE,MAGENTA,CYAN,VIOLET,PINK,ORANGE,
                YELLOW,SPRING,TEAL,AZURE,DARKRED,DARKGREEN,DARKBLUE,DARKMAGENTA,
                DARKCYAN,DARKVIOLET,DARKPINK,DARKORANGE,DARKYELLOW,LIGHTGREEN,
                LIGHTBLUE,LIGHTRED,LIGHTMAGENTA,LIGHTCYAN,LIGHTVIOLET,LIGHTPINK,
                LIGHTORANGE,LIGHTYELLOW};

            for (const auto mc : mcPrimaryVector)
            {
                if(primaryHits.find(mc) != primaryHits.end())
                {
                    const std::string primaryName(viewStr + "PrimaryHits_" + std::to_string(iPrimary));

                    if (iPrimary < 31)
                    {
                        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                    &primaryHits[mc], primaryName, colors[iPrimary]));
                    }
                    else
                    {
                        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                    &primaryHits[mc], primaryName, WHITE));
                    }
                }
                iPrimary++;
            }
        }
    }

    if (m_tagPid || m_tagPrimary)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------

StatusCode VisualHitTaggingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "TagPID", m_tagPid));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "TagPrimary", m_tagPrimary));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "RetainNonReconstructable", m_retainNonReconstructable));

    return STATUS_CODE_SUCCESS;
}

/*
StatusCode VisualHitTaggingAlgorithm::Train()
{
    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW) return STATUS_CODE_NOT_ALLOWED;

        std::string trainingOutputFileName(m_trainingOutputFile);
        LArMvaHelper::MvaFeatureVector featureVector;

        if (isU) trainingOutputFileName += "_CaloHitListU.txt";
        else if (isV) trainingOutputFileName += "_CaloHitListV.txt";
        else if (isW) trainingOutputFileName += "_CaloHitListW.txt";

        featureVector.push_back(static_cast<double>(pCaloHitList->size()));

        int discardedHits = 0;
        LArMCParticleHelper::PrimaryParameters parameters;
        // Only care about reconstructability with respect to the current view, so skip good view check
        parameters.m_minHitsForGoodView = 0;
        // Turn off max photo propagation for now, only care about killing off daughters of neutrons
        parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList,
                parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            int isReconstructable = 1;
            int pdg(-1);

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                // Throw away non-reconstructable hits
                if (targetMCParticleToHitsMap.find(pMCParticle) == targetMCParticleToHitsMap.end())
                {
                    ++discardedHits;
                    isReconstructable = 0;
                }
                pdg = pMCParticle->GetParticleId();
                if (isReconstructable)
                {
                    if(isNeutronDaughter(pMCParticle))
                    {
                        isReconstructable = 0;
                    }
                }
            }
            catch (...)
            {
                continue;
            }

            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetY()));
            featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
            featureVector.push_back(static_cast<double>(pdg));
            featureVector.push_back(static_cast<double>(isReconstructable));
        }

        LArMvaHelper::ProduceTrainingExample(trainingOutputFileName, true, featureVector);
    }
    return STATUS_CODE_SUCCESS;
}
*/

} // namespace lar_content
