/**
 *  @file   larpandoracontent/LArTwoDReco/CheatingStreamSelectionAlgorithm.cc
 *
 *  @brief  Implementation of cheating stream selection algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingStreamSelectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

StatusCode CheatingStreamSelectionAlgorithm::AllocateToStreams(const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList{pCluster->GetOrderedCaloHitList()};
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);
    const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
    caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
    FloatVector trackAdc, showerAdc;
    for (const CaloHit *pCaloHit : caloHits)
    {
        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const float adc{pCaloHit->GetInputEnergy()};
            if (adc < 0.f)
                continue;

            const int pdg{std::abs(pMCParticle->GetParticleId())};
            if (pdg == 11 || pdg == 22)
                showerAdc.emplace_back(adc);
            else
                trackAdc.emplace_back(adc);
        }
        catch (const StatusCodeException &)
        {
            continue;
        }
    }

    const float trackAcc{std::accumulate(std::begin(trackAdc), std::end(trackAdc), 0.f)};
    const float showerAcc{std::accumulate(std::begin(showerAdc), std::end(showerAdc), 0.f)};

    if (trackAcc >= showerAcc)
        m_clusterListMap.at(m_trackListName).emplace_back(pCluster);
    else
        m_clusterListMap.at(m_showerListName).emplace_back(pCluster);

    return STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingStreamSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, StreamSelectionAlgorithm::ReadSettings(xmlHandle));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackListName", m_trackListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerListName", m_showerListName));

    m_listNames.emplace_back(m_trackListName);
    m_listNames.emplace_back(m_showerListName);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
