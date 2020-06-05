/**
 *  @file   larpandoracontent/LArTrackShowerId/NnClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTrackShowerId/NnClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

bool NnClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() <= 25) // Network not reliable for low hits clusters, keep existing Id
        return pCluster->GetParticleId() == MU_MINUS;

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    FloatVector probabilities;
    float accum{0.f};
    for (const CaloHit *pCaloHit : caloHitList)
    {
        auto properties = pCaloHit->GetPropertiesMap();
        const float pTrack{properties["Ptrack"]};
        probabilities.push_back(pTrack);
        accum += pTrack;
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
            return pCluster->GetParticleId() == MU_MINUS;
    }
    else
    {
        const float critical{mu + sigma};
        if (critical < 0.5f)    // Track
            return false;
        else                    // Ambiguous, go with the current Id
            return pCluster->GetParticleId() == MU_MINUS;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NnClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ClusterCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
