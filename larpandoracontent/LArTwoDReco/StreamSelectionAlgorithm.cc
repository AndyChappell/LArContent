/**
 *  @file   larpandoradlcontent/LArTwoDReco/StreamSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/StreamSelectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

StreamSelectionAlgorithm::StreamSelectionAlgorithm() :
    m_useTracksAsOutputList{true}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StreamSelectionAlgorithm::~StreamSelectionAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StreamSelectionAlgorithm::Run()
{
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    std::string originalClusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalClusterListName));

    ClusterList trackClusterList, showerClusterList;
    for (const Cluster *pCluster : *pClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList{pCluster->GetOrderedCaloHitList()};
        CaloHitList caloHits;
        orderedCaloHitList.FillCaloHitList(caloHits);
        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
        FloatVector trackLikelihoods;
        try
        {
            for (const CaloHit *pCaloHit : caloHits)
            {
                const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
                const float pTrack{pLArCaloHit->GetTrackProbability()};
                const float pShower{pLArCaloHit->GetShowerProbability()};
                if ((pTrack + pShower) > std::numeric_limits<float>::epsilon())
                    trackLikelihoods.emplace_back(pTrack / (pTrack + pShower));
            }

            const unsigned long N{trackLikelihoods.size()};
            if (N > 0)
            {
                float mean{std::accumulate(std::begin(trackLikelihoods), std::end(trackLikelihoods), 0.f) / N};
                /*float accum{0.0};
                std::for_each (std::begin(trackLikelihoods), std::end(trackLikelihoods), [&](const float x)
                    {
                        accum += (x - mean) * (x - mean);
                    });

                float stdev{N > 1 ? std::sqrt(accum / (N - 1)) : 0.f};*/
                if (mean >= 0.5f)
                    trackClusterList.emplace_back(pCluster);
                else
                    showerClusterList.emplace_back(pCluster);
            }
        }
        catch (const StatusCodeException&)
        {
            continue;
        }
    }

    // ATTN - We're ok with saving empty lists here and allowing future algorithms to simply do nothing if there are no clusters
    // Moves the subset of clusters in the cluster list from the old list to the new list
    PandoraContentApi::SaveList<ClusterList>(*this, originalClusterListName, m_trackClusterListName, trackClusterList);
    PandoraContentApi::SaveList<ClusterList>(*this, originalClusterListName, m_showerClusterListName, showerClusterList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StreamSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackClusterListName", m_trackClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerClusterListName", m_showerClusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TracksAreOutput",
        m_useTracksAsOutputList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
