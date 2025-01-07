/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

KalmanClusterCreationAlgorithm::KalmanClusterCreationAlgorithm() :
    m_minMipFraction{0.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::Run()
{
    for (std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (!pCaloHitList->empty())
        {
            HitType view{pCaloHitList->front()->GetHitType()};
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FilterCaloHits(pCaloHitList, m_viewHitsMap[view]));
        }
    }

    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        size_t totalHits{0};
        for (OrderedCaloHitList::const_iterator iter = m_viewHitsMap[view].begin(), iterEnd = m_viewHitsMap[view].end(); iter != iterEnd; ++iter)
            totalHits += iter->second->size();
        std::cout << "View: " << view << " Filtered: " << totalHits << std::endl;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::FilterCaloHits(const CaloHitList *const pCaloHitList, OrderedCaloHitList &selectedCaloHitList) const
{
    CaloHitList availableHitList;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (PandoraContentApi::IsAvailable(*this, pCaloHit) && pCaloHit->GetMipEquivalentEnergy() >= m_minMipFraction)
            availableHitList.push_back(pCaloHit);
    }

    if (availableHitList.empty())
        return STATUS_CODE_SUCCESS;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, selectedCaloHitList.Add(availableHitList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMipFraction", m_minMipFraction));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
