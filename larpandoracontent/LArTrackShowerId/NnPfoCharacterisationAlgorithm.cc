/**
 *  @file   larpandoracontent/LArTrackShowerId/NnPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArTrackShowerId/NnPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

bool NnPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
    FloatVector probabilities;
    float accum{0.f};
    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (const Cluster *pCluster : pfoClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
        for (const CaloHit *pCaloHit : caloHitList)
        {
            auto properties = pCaloHit->GetPropertiesMap();
            const float pTrack{properties["Ptrack"]};
            probabilities.push_back(pTrack);
            accum += pTrack;
        }
    }

    const size_t N{probabilities.size()};
    const float mu{accum / probabilities.size()};

    accum = 0.f;
    std::for_each(probabilities.begin(), probabilities.end(), [&](const float p){
        const float diff{p - mu};
        accum += diff * diff;
    });
    const float sigma{std::sqrt(accum / (N - 1))};
    if (mu > 0.5f)
    {
        const float critical{mu - sigma};
        if (critical > 0.5f)    // Track
            return true;
        else                    // Ambiguous, go with the current Id
            return pPfo->GetParticleId() == MU_MINUS;
    }
    else
    {
        const float critical{mu + sigma};
        if (critical < 0.5f)    // Track
            return false;
        else                    // Ambiguous, go with the current Id
            return pPfo->GetParticleId() == MU_MINUS;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NnPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const /*pCluster*/) const
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NnPfoCharacterisationAlgorithm::IsClearTrack3x2D(const ParticleFlowObject *const /*pPfo*/) const
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NnPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
