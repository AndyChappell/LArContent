/**
 *  @file   larpandoradlcontent/LArMonitoring/DlClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower cluster validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArMonitoring/DlClusterValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlClusterValidationAlgorithm::DlClusterValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlClusterValidationAlgorithm::~DlClusterValidationAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_outputTreeName, m_outputFileName, "UPDATE"));
    }
    catch(const StatusCodeException&)
    {
        std::cout << "DlClusterValidationAlgorithm: Unable to write tree to file" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterValidationAlgorithm::Run()
{
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 0;
    parameters.m_minHitsForGoodView = 0;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_minHitSharingFraction = 0.f;
    parameters.m_foldBackHierarchy = false;
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    for (const Cluster *pCluster : *pClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList{pCluster->GetOrderedCaloHitList()};
        CaloHitList caloHits;
        orderedCaloHitList.FillCaloHitList(caloHits);
        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
        FloatVector trackLikelihoods;
        std::map<const MCParticle*, int> mcHits;
        try
        {
            for (const CaloHit *pCaloHit : caloHits)
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                if (!pMCParticle)
                    continue;
                if (mcHits.find(pMCParticle) != mcHits.end())
                    ++mcHits[pMCParticle];
                else
                    mcHits[pMCParticle] = 1;

                const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
                const float pTrack{pLArCaloHit->GetTrackProbability()};
                const float pShower{pLArCaloHit->GetShowerProbability()};
                if ((pTrack + pShower) > std::numeric_limits<float>::epsilon())
                    trackLikelihoods.emplace_back(pTrack / (pTrack + pShower));
            }

            // Find best MC truth match
            const MCParticle *bestMC{nullptr};
            int bestCount{0};
            for (const auto [ key, value ] : mcHits)
            {
                if (value > bestCount)
                {
                    bestMC = key;
                    bestCount = value;
                }
            }
            if (!bestMC)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);
            const int pdg{std::abs(bestMC->GetParticleId())};

            const unsigned long N{trackLikelihoods.size()};
            if (N > 1)
            {
                float mean{std::accumulate(std::begin(trackLikelihoods), std::end(trackLikelihoods), 0.f) / N};
                /*float accum{0.0};
                std::for_each (std::begin(trackLikelihoods), std::end(trackLikelihoods), [&](const float x)
                    {
                        accum += (x - mean) * (x - mean);
                    });

                float stdev{std::sqrt(accum / (N - 1))};*/
                int trueTrack{0};
                int falseTrack{0};
                int trueShower{0};
                int falseShower{0};
                if (mean >= 0.5f)
                {   // Track classification
                    if (pdg == 11 || pdg == 22)
                        ++falseTrack;
                    else
                        ++trueTrack;
                }
                else
                {   // Shower classification
                    if (pdg == 11 || pdg == 22)
                        ++trueShower;
                    else
                        ++falseShower;
                }
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "true_track", trueTrack));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "false_track", falseTrack));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "true_shower", trueShower));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "false_shower", falseShower));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_outputTreeName));
            }
        }
        catch (const StatusCodeException&)
        {
            continue;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFileName", m_outputFileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTreeName", m_outputTreeName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
