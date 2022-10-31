/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/HitCorrelationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCorrelationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

HitCorrelationAlgorithm::HitCorrelationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitCorrelationAlgorithm::~HitCorrelationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCorrelationAlgorithm::Run()
{
    m_volumeMap.clear();

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList, m_caloHitListName));

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        if (!pLArCaloHit)
            continue;
        const unsigned int tpcVolume{pLArCaloHit->GetLArTPCVolumeId()};
        const unsigned int childVolume{pLArCaloHit->GetDaughterVolumeId()};
        const unsigned int key{1000 * tpcVolume + childVolume};
        m_volumeMap[key].AddCaloHit(pLArCaloHit);
        m_availabilityMap[pLArCaloHit] = true;
    }

    for (const auto & [ key, volume ] : m_volumeMap)
    {
        std::cout << "Volume " << key << std::endl; 
        const CaloHitList caloHitListU{volume.GetCaloHits(HitType::TPC_VIEW_U)};
        const CaloHitList caloHitListV{volume.GetCaloHits(HitType::TPC_VIEW_V)};
        const CaloHitList caloHitListW{volume.GetCaloHits(HitType::TPC_VIEW_W)};
        HitMap hitMap;
        this->Correlate(caloHitListU, caloHitListV, hitMap);
        this->Correlate(caloHitListU, caloHitListW, hitMap);
        this->Correlate(caloHitListV, caloHitListW, hitMap);

        HitTable usedHits;
        LArTripletVector hitTriplets;
        for (const auto & [pCaloHit, caloHits] : hitMap)
        {
            LArSet hitSet{pCaloHit};
            this->FindRelationships(caloHits, hitMap, hitSet);
            this->MakeHitTriplets(hitSet, usedHits, hitTriplets);
        }

        for (const LArTriplet &triplet : hitTriplets)
        {
            const CaloHit *const pCaloHitU{std::get<0>(triplet)};
            const CaloHit *const pCaloHitV{std::get<1>(triplet)};
            const CaloHit *const pCaloHitW{std::get<2>(triplet)};
            const CartesianVector &posU{pCaloHitU->GetPositionVector()};
            const CartesianVector &posV{pCaloHitV->GetPositionVector()};
            const CartesianVector &posW{pCaloHitW->GetPositionVector()};
            CartesianVector pos3D(0, 0, 0);
            float chi2{0.f};
            LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, posU, posV, posW, pos3D, chi2);
            std::cout << "In [(" << posU.GetX() << "," << posU.GetZ() << ") (" << posV.GetX() << "," << posV.GetZ() << ") (" <<
                posW.GetX() << "," << posW.GetZ() << ")] Out (" << pos3D.GetX() << "," << pos3D.GetZ() << ") chi2: " << chi2 << std::endl;
        }
        std::cout << "Inputs: " << caloHitListU.size() << " " << caloHitListV.size() << " " << caloHitListW.size() << std::endl;
        std::cout << "Triplets: " << hitTriplets.size() << std::endl;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCorrelationAlgorithm::FindRelationships(const CaloHitList &caloHitList, const HitMap &hitMap, LArSet &hitSet) const
{
    for (const CaloHit *pCaloHit : caloHitList)
    {
        if (hitSet.find(pCaloHit) == hitSet.end())
        {
            hitSet.insert(pCaloHit);
            this->FindRelationships(hitMap.at(pCaloHit), hitMap, hitSet);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCorrelationAlgorithm::Correlate(const CaloHitList &caloHitList1, const CaloHitList &caloHitList2, HitMap &hitMap) const
{
    if (caloHitList1.empty() || caloHitList2.empty())
        return;

    const float minX{std::max(caloHitList1.front()->GetPositionVector().GetX(), caloHitList2.front()->GetPositionVector().GetX())};
    const float maxX{std::min(caloHitList1.back()->GetPositionVector().GetX(), caloHitList2.back()->GetPositionVector().GetX())};

    if (minX > maxX)
        return;

    std::cout << "Overlap found in " << caloHitList1.front()->GetHitType() << " and " << caloHitList2.front()->GetHitType() << " from " <<
        minX << " to " << maxX << std::endl;
    auto iter1{caloHitList1.begin()};
    while((*iter1)->GetPositionVector().GetX() < minX)
        ++iter1;

    while (iter1 != caloHitList1.end())
    {
        const CaloHit *const pCaloHit1{*iter1};
        const float x1{pCaloHit1->GetPositionVector().GetX()};
        if (x1 > maxX)
            break;
        // Make this half the wire pitch for W
        const float lowX{x1 - 0.025f};
        const float highX{x1 + 0.025f};

        for (const CaloHit *pCaloHit2 : caloHitList2)
        {
            const float x2{pCaloHit2->GetPositionVector().GetX()};
            if ((x2 >= lowX) && (x2 <= highX))
            {
                hitMap[pCaloHit1].emplace_back(pCaloHit2);
                hitMap[pCaloHit2].emplace_back(pCaloHit1);
            }
        }
        ++iter1;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCorrelationAlgorithm::MakeHitTriplets(const LArSet &caloHitSet, HitTable &usedHits, LArTripletVector &hitTriplets) const
{
    CaloHitList caloHitListU, caloHitListV, caloHitListW;

    for (const CaloHit *const pCaloHit : caloHitSet)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                caloHitListU.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_V:
                caloHitListV.emplace_back(pCaloHit);
                break;
            default:
                caloHitListW.emplace_back(pCaloHit);
                break;
        }
    }

    CartesianVector pos3D(0, 0, 0);
    for (const CaloHit *const pCaloHitU : caloHitListU)
    {
        if (usedHits.find(pCaloHitU) != usedHits.end())
            continue;
        LArTriplet bestTriplet({nullptr, nullptr, nullptr});
        float bestChi2{std::numeric_limits<float>::max()};
        const CartesianVector &posU{pCaloHitU->GetPositionVector()};
        for (const CaloHit *const pCaloHitV : caloHitListV)
        {
            if (usedHits.find(pCaloHitV) != usedHits.end())
                continue;
            const CartesianVector &posV{pCaloHitV->GetPositionVector()};
            for (const CaloHit *const pCaloHitW : caloHitListW)
            {
                if (usedHits.find(pCaloHitW) != usedHits.end())
                    continue;
                const CartesianVector &posW{pCaloHitW->GetPositionVector()};
                float chi2{0.f};
                LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, posU, posV, posW, pos3D, chi2);
                if (chi2 < 0.1f)
                {
                    if (chi2 < bestChi2)
                    {
                        bestTriplet = std::make_tuple(pCaloHitU, pCaloHitV, pCaloHitW);
                        bestChi2 = chi2;
                    }
                }
            }
        }
        if (bestChi2 < std::numeric_limits<float>::max())
        {
            usedHits[std::get<0>(bestTriplet)] = true;
            usedHits[std::get<1>(bestTriplet)] = true;
            usedHits[std::get<2>(bestTriplet)] = true;
            hitTriplets.emplace_back(bestTriplet);
        }
    }

    // Need to consider leftover hits and whether they can form 2 hit combos
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCorrelationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void HitCorrelationAlgorithm::TpcChildVolume::AddCaloHit(const CaloHit *pCaloHit)
{
    switch (pCaloHit->GetHitType())
    {
        case HitType::TPC_VIEW_U:
            m_hitsU.emplace_back(pCaloHit);
            break;
        case HitType::TPC_VIEW_V:
            m_hitsV.emplace_back(pCaloHit);
            break;
        default:
            m_hitsW.emplace_back(pCaloHit);
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList HitCorrelationAlgorithm::TpcChildVolume::GetCaloHits(const HitType type) const
{
    switch (type)
    {
        case HitType::TPC_VIEW_U:
        {
            std::vector<const CaloHit *> sortedHitList;
            for (const CaloHit *pCaloHit : m_hitsU)
                sortedHitList.emplace_back(pCaloHit);
            std::sort(sortedHitList.begin(), sortedHitList.end(), LArClusterHelper::SortHitsByPositionInX);
            CaloHitList caloHitList;
            for (const CaloHit *pCaloHit : sortedHitList)
                caloHitList.emplace_back(pCaloHit);
            return caloHitList;
        }
        case HitType::TPC_VIEW_V:
        {
            std::vector<const CaloHit *> sortedHitList;
            for (const CaloHit *pCaloHit : m_hitsV)
                sortedHitList.emplace_back(pCaloHit);
            std::sort(sortedHitList.begin(), sortedHitList.end(), LArClusterHelper::SortHitsByPositionInX);
            CaloHitList caloHitList;
            for (const CaloHit *pCaloHit : sortedHitList)
                caloHitList.emplace_back(pCaloHit);
            return caloHitList;
        }
        default:
        {
            std::vector<const CaloHit *> sortedHitList;
            for (const CaloHit *pCaloHit : m_hitsW)
                sortedHitList.emplace_back(pCaloHit);
            std::sort(sortedHitList.begin(), sortedHitList.end(), LArClusterHelper::SortHitsByPositionInX);
            CaloHitList caloHitList;
            for (const CaloHit *pCaloHit : sortedHitList)
                caloHitList.emplace_back(pCaloHit);
            return caloHitList;
        }
    }
}

} // namespace lar_content
