/**
 *  @file   larpandoracontent/LArMonitoring/PiZeroMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/PiZeroMonitoringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PiZeroMonitoringAlgorithm::PiZeroMonitoringAlgorithm() :
    m_caloHitListName("CaloHitList2D"),
    m_pfoListName("RecreatedPfos")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PiZeroMonitoringAlgorithm::~PiZeroMonitoringAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "pizero", "pizero.root", "UPDATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "PiZeroMonitoringAlgorithm: Unable to write to ROOT tree" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PiZeroMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    LArMCParticleHelper::MCContributionMap mcHitMap;
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            mcHitMap[pMC].emplace_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
            continue;
        }
    }

    const Pfo* pLongestTrackPfo{nullptr};
    float longestTrackLength{0.f};
    std::unordered_map<const Pfo *, const MCParticle *> pfoMCMap;
    std::unordered_map<const Pfo*, CaloHitList> pfoHitMap;
    CartesianVector nuVertex{0.f, 0.f, 0.f};
    for (const Pfo* pPfo : *pPfoList)
    {
        const int pfoPdg{std::abs(pPfo->GetParticleId())};
        if (pfoPdg == 12 || pfoPdg == 14 || pfoPdg == 16)
        {
            const Vertex *pVtx{LArPfoHelper::GetVertex(pPfo)};
            nuVertex.SetValues(pVtx->GetPosition().GetX(), pVtx->GetPosition().GetY(), pVtx->GetPosition().GetZ());
            continue;
        }

        ClusterList cluster3DList;
        LArPfoHelper::GetClusters(pPfo, TPC_3D, cluster3DList);
        const Cluster *pCluster3D{cluster3DList.empty() ? nullptr : cluster3DList.front()};
        // Only consider PFOs with 3D clusters
        if (pCluster3D == nullptr)
            continue;

        std::unordered_map<const MCParticle*, float> mcAdcMap;
        CaloHitList all2DHits;
        LArPfoHelper::GetAllCaloHits(pPfo, all2DHits);
        for (const CaloHit* pCaloHit : all2DHits)
        {
            try
            {
                const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (mcAdcMap.find(pMC) == mcAdcMap.end())
                    mcAdcMap[pMC] = 0.f;
                mcAdcMap[pMC] += pCaloHit->GetInputEnergy();
                pfoHitMap[pPfo].emplace_back(pCaloHit);
            }
            catch (StatusCodeException &)
            {
                continue;
            }
        }
        const MCParticle *pBestMatchMC{nullptr};
        float bestMatchAdc{0.f};
        for (const auto& [pMC, adc] : mcAdcMap)
        {
            if (adc > bestMatchAdc)
            {
                bestMatchAdc = adc;
                pBestMatchMC = pMC;
            }
        }
        pfoMCMap[pPfo] = pBestMatchMC;

        const float length{LArClusterHelper::GetLength(pCluster3D)};
        if (length > longestTrackLength)
        {
            longestTrackLength = length;
            pLongestTrackPfo = pPfo;
        }
    }

    for (const auto& [pPfo, pMC] : pfoMCMap)
    {
        CaloHitList pfoHitsU, pfoHitsV, pfoHitsW;
        for (const CaloHit* pCaloHit : pfoHitMap[pPfo])
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    pfoHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    pfoHitsV.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    pfoHitsW.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
        CaloHitList mcHitsU, mcHitsV, mcHitsW;
        for (const CaloHit* pCaloHit : mcHitMap[pMC])
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    mcHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    mcHitsV.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    mcHitsW.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
/*
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
        if (!mcHitsU.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mcHitsU, "MC U", RED));
        }
        if (!pfoHitsU.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoHitsU, "PFO U", ORANGE));
        }
        if (!mcHitsV.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mcHitsV, "MC V", VIOLET));
        }
        if (!pfoHitsV.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoHitsV, "PFO V", PINK));
        }
        if (!mcHitsW.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mcHitsW, "MC W", BLUE));
        }
        if (!pfoHitsW.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &pfoHitsW, "PFO W", CYAN));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
*/
        std::cout << "PFO " << pPfo << " matched to MC particle " << pMC << " with PDG " << pMC->GetParticleId() << std::endl;
    }

    // Get the ADC counts for each shower candidate
    std::map<float, const Pfo *> adcPfoMap;
    for (const auto& [pPfo, pMC] : pfoMCMap)
    {
        if (pPfo != pLongestTrackPfo)
        {
            float adcSumU{0.f}, adcSumV{0.f}, adcSumW{0.f};
            for (const CaloHit *const pCaloHit : pfoHitMap[pPfo])
            {
                switch (pCaloHit->GetHitType())
                {
                    case TPC_VIEW_U:
                        adcSumU += pCaloHit->GetInputEnergy();
                        break;
                    case TPC_VIEW_V:
                        adcSumV += pCaloHit->GetInputEnergy();
                        break;
                    case TPC_VIEW_W:
                        adcSumW += pCaloHit->GetInputEnergy();
                        break;
                    default:
                        break;
                }
            }
            adcPfoMap[std::max({adcSumU, adcSumV, adcSumW})] = pPfo;
        }
    }
    if (adcPfoMap.size() >= 2)
    {
        // Pick the twos PFOs with the highest ADC counts as the shower candidates
        const Pfo *pShowerCandidate1{adcPfoMap.rbegin()->second};
        const Pfo *pShowerCandidate2{std::next(adcPfoMap.rbegin())->second};
        std::cout << "Shower candidate 1: PFO " << pShowerCandidate1 << " with ADC " << adcPfoMap.rbegin()->first << std::endl;
        std::cout << "Shower candidate 2: PFO " << pShowerCandidate2 << " with ADC " << std::next(adcPfoMap.rbegin())->first << std::endl;

        auto cmp = [](const CaloHit* a, const CaloHit* b) {
            return a < b;
        };
        CaloHitList pfo1Hits(pfoHitMap[pShowerCandidate1]);
        CaloHitList mc1Hits(mcHitMap[pfoMCMap[pShowerCandidate1]]);
        pfo1Hits.sort(cmp);
        mc1Hits.sort(cmp);
        CaloHitList pfo2Hits(pfoHitMap[pShowerCandidate2]);
        CaloHitList mc2Hits(mcHitMap[pfoMCMap[pShowerCandidate2]]);
        pfo2Hits.sort(cmp);
        mc2Hits.sort(cmp);

        CaloHitList intersection1, intersection2;
        std::set_intersection(mc1Hits.begin(), mc1Hits.end(), pfo1Hits.begin(), pfo1Hits.end(), std::back_inserter(intersection1));
        std::set_intersection(mc2Hits.begin(), mc2Hits.end(), pfo2Hits.begin(), pfo2Hits.end(), std::back_inserter(intersection2));

        CaloHitList shower13DHits1, shower13DHits2;
        LArPfoHelper::GetCaloHits(pShowerCandidate1, TPC_3D, shower13DHits1);
        LArPfoHelper::GetCaloHits(pShowerCandidate2, TPC_3D, shower13DHits2);
        CartesianVector shower1Dir{0.f, 0.f, 0.f}, shower2Dir{0.f, 0.f, 0.f};
        float weight1{0.f}, weight2{0.f};
        for (const CaloHit* pCaloHit : shower13DHits1)
        {
            const CartesianVector localDir{(pCaloHit->GetPositionVector() - nuVertex) * pCaloHit->GetInputEnergy()};
            shower1Dir += localDir;
            weight1 += pCaloHit->GetInputEnergy();
        }
        for (const CaloHit* pCaloHit : shower13DHits2)
        {
            const CartesianVector localDir{(pCaloHit->GetPositionVector() - nuVertex) * pCaloHit->GetInputEnergy()};
            shower2Dir += localDir;
            weight2 += pCaloHit->GetInputEnergy();
        }
        shower1Dir *= 1.f / weight1;
        shower2Dir *= 1.f / weight2;
        shower1Dir = shower1Dir.GetUnitVector();
        shower2Dir = shower2Dir.GetUnitVector();

        const MCParticle *pMC1{pfoMCMap[pShowerCandidate1]};
        const MCParticle *pMC2{pfoMCMap[pShowerCandidate2]};
        const CartesianVector mc1Dir{pMC1->GetMomentum().GetUnitVector()};
        const CartesianVector mc2Dir{pMC2->GetMomentum().GetUnitVector()};
        std::cout << "Reco dir 1: " << shower1Dir.GetX() << ", " << shower1Dir.GetY() << ", " << shower1Dir.GetZ() << std::endl;
        std::cout << "MC dir 1: " << mc1Dir.GetX() << ", " << mc1Dir.GetY() << ", " << mc1Dir.GetZ() << std::endl;
        std::cout << "Reco dir 2: " << shower2Dir.GetX() << ", " << shower2Dir.GetY() << ", " << shower2Dir.GetZ() << std::endl;
        std::cout << "MC dir 2: " << mc2Dir.GetX() << ", " << mc2Dir.GetY() << ", " << mc2Dir.GetZ() << std::endl;

        const float purity1{static_cast<float>(intersection1.size()) / pfoHitMap[pShowerCandidate1].size()};
        const float completeness1{static_cast<float>(intersection1.size()) / mcHitMap[pfoMCMap[pShowerCandidate1]].size()};
        const float purity2{static_cast<float>(intersection2.size()) / pfoHitMap[pShowerCandidate2].size()};
        const float completeness2{static_cast<float>(intersection2.size()) / mcHitMap[pfoMCMap[pShowerCandidate2]].size()};
        std::cout << "Shower candidate 1: purity = " << purity1 << ", completeness = " << completeness1 << std::endl;
        std::cout << "Shower candidate 2: purity = " << purity2 << ", completeness = " << completeness2 << std::endl;

        float adcFrac1{0.f}, adcFrac2{0.f}, mcAdcFrac1{0.f}, mcAdcFrac2{0.f};
        for (const CaloHit* pCaloHit : pfo1Hits)
            adcFrac1 += pCaloHit->GetInputEnergy();
        for (const CaloHit* pCaloHit : pfo2Hits)
            adcFrac2 += pCaloHit->GetInputEnergy();
        for (const CaloHit* pCaloHit : mc1Hits)
            mcAdcFrac1 += pCaloHit->GetInputEnergy();
        for (const CaloHit* pCaloHit : mc2Hits)
            mcAdcFrac2 += pCaloHit->GetInputEnergy();
        adcFrac1 /= mcAdcFrac1;
        adcFrac2 /= mcAdcFrac2;
        const int efficient{pMC1 != pMC2};
        const float openingAngle1{static_cast<float>(180 / M_PI) * std::acos(shower1Dir.GetDotProduct(mc1Dir))};
        const float openingAngle2{static_cast<float>(180 / M_PI) * std::acos(shower2Dir.GetDotProduct(mc2Dir))};

        const float invariantMassTrueE{std::sqrt(2 * pMC1->GetEnergy() * pMC2->GetEnergy() * (1 - shower1Dir.GetDotProduct(shower2Dir)))};
        const float invariantMassRecoE{std::sqrt(2 * adcFrac1 * pMC1->GetEnergy() * adcFrac2 * pMC2->GetEnergy() * (1 - shower1Dir.GetDotProduct(shower2Dir)))};
        std::cout << "Invariant mass: " << invariantMassTrueE << " GeV (true energy), " << invariantMassRecoE << " GeV (reco energy)" << std::endl;

        // Populate ROOT tree
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "efficient", efficient));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "openingAngle1", openingAngle1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "openingAngle2", openingAngle2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "purity1", purity1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "completeness1", completeness1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "purity2", purity2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "completeness2", completeness2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "adcFrac1", adcFrac1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "adcFrac2", adcFrac2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "invariantMassTrueE", invariantMassTrueE));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "invariantMassRecoE", invariantMassRecoE));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "pizero"));
    }
    else if (adcPfoMap.size() == 1)
    {
        const Pfo *pShowerCandidate1{adcPfoMap.rbegin()->second};
        std::cout << "Shower candidate 1: PFO " << pShowerCandidate1 << " with ADC " << adcPfoMap.rbegin()->first << std::endl;

        auto cmp = [](const CaloHit* a, const CaloHit* b) {
            return a < b;
        };

        CaloHitList pfo1Hits(pfoHitMap[pShowerCandidate1]);
        CaloHitList mc1Hits(mcHitMap[pfoMCMap[pShowerCandidate1]]);
        pfo1Hits.sort(cmp);
        mc1Hits.sort(cmp);

        const MCParticle *pMC1{pfoMCMap[pShowerCandidate1]};
        const CartesianVector mc1Dir{pMC1->GetMomentum().GetUnitVector()};

        CaloHitList shower13DHits1;
        LArPfoHelper::GetCaloHits(pShowerCandidate1, TPC_3D, shower13DHits1);
        CartesianVector shower1Dir{0.f, 0.f, 0.f};
        float weight1{0.f};
        for (const CaloHit* pCaloHit : shower13DHits1)
        {
            const CartesianVector localDir{(pCaloHit->GetPositionVector() - nuVertex) * pCaloHit->GetInputEnergy()};
            shower1Dir += localDir;
            weight1 += pCaloHit->GetInputEnergy();
        }
        shower1Dir *= 1.f / weight1;
        shower1Dir = shower1Dir.GetUnitVector();

        const float openingAngle1{static_cast<float>(180 / M_PI) * std::acos(shower1Dir.GetDotProduct(mc1Dir))};
        const float openingAngle2{-1};

        float adcFrac1{0.f}, mcAdcFrac1{0.f};
        for (const CaloHit* pCaloHit : pfo1Hits)
            adcFrac1 += pCaloHit->GetInputEnergy();
        for (const CaloHit* pCaloHit : mc1Hits)
            mcAdcFrac1 += pCaloHit->GetInputEnergy();
        adcFrac1 /= mcAdcFrac1;

        CaloHitList intersection1;
        std::set_intersection(mc1Hits.begin(), mc1Hits.end(), pfo1Hits.begin(), pfo1Hits.end(), std::back_inserter(intersection1));
        const float purity1{static_cast<float>(intersection1.size()) / pfoHitMap[pShowerCandidate1].size()};
        const float completeness1{static_cast<float>(intersection1.size()) / mcHitMap[pfoMCMap[pShowerCandidate1]].size()};
        const float purity2{-1}, completeness2{-1}, adcFrac2{-1}, invariantMassTrueE{-1}, invariantMassRecoE{-1};
        const int efficient{0};

        std::cout << "Shower candidate 1: purity = " << purity1 << ", completeness = " << completeness1 << std::endl;

        std::cout << "Not enough shower candidates (" << adcPfoMap.size() << ") found." << std::endl;

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "efficient", efficient));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "openingAngle1", openingAngle1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "openingAngle2", openingAngle2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "purity1", purity1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "completeness1", completeness1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "purity2", purity2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "completeness2", completeness2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "adcFrac1", adcFrac1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "adcFrac2", adcFrac2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "invariantMassTrueE", invariantMassTrueE));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "invariantMassRecoE", invariantMassRecoE));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "pizero"));
    }
    else
    {
        const int efficient{0};
        const float openingAngle1{-1}, openingAngle2{-1}, purity1{-1}, completeness1{-1}, purity2{-1}, completeness2{-1}, adcFrac1{-1}, adcFrac2{-1}, invariantMassTrueE{-1}, invariantMassRecoE{-1};
        std::cout << "No shower candidates found." << std::endl;
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "efficient", efficient));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "openingAngle1", openingAngle1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "openingAngle2", openingAngle2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "purity1", purity1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "completeness1", completeness1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "purity2", purity2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "completeness2", completeness2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "adcFrac1", adcFrac1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "adcFrac2", adcFrac2));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "invariantMassTrueE", invariantMassTrueE));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "pizero", "invariantMassRecoE", invariantMassRecoE));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "pizero"));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PiZeroMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
