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
        const CaloHitList caloHitListU{volume.GetCaloHits(HitType::TPC_VIEW_U)};
        const CaloHitList caloHitListV{volume.GetCaloHits(HitType::TPC_VIEW_V)};
        const CaloHitList caloHitListW{volume.GetCaloHits(HitType::TPC_VIEW_W)};
        HitMap hitMap;
        this->Correlate(caloHitListU, caloHitListV, hitMap);
        this->Correlate(caloHitListU, caloHitListW, hitMap);
        this->Correlate(caloHitListV, caloHitListW, hitMap);

        // Having found the 2D pairs, find the links to the third view and create sets of hits that ought to be considered together
    }

    return STATUS_CODE_SUCCESS;
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
    // Next we need to figure out an efficient way to step through the lists and group hits (ideally we'd like to note where we are
    // in the list and step back until out of range
    //
    // e.g. Hit 1 list 1 has a position, find the first hit in the other list within bounds and walk forward until out of bounds
    // then Hit 2 list 1 has a position, walk back (if necessary) from where we were in list 2 until we hit the lower limit, repeat
    // This should guarantee that all hits in list 2 should also have their corresponding partners
    (void)hitMap;
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
