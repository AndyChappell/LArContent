/**
 *  @file   larpandoracontent/LArTrackShowerId/NnPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include <cmath>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArTrackShowerId/NnPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NnPfoCharacterisationAlgorithm::NnPfoCharacterisationAlgorithm() 
{
}

//-----------------------------------------------------------------------------
float GetTrackLikelihood(const Cluster *const pCluster)
{
    const OrderedCaloHitList& orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float trackProb{0.f}, showerProb{0.f};
    for(const auto caloHitList : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *caloHitList.second)
        {
            auto properties = pCaloHit->GetPropertiesMap();
            trackProb += std::log(properties["Ptrack"] + 1e-5);
            showerProb += std::log(properties["Pshower"] + 1e-5);
        }
    }
    trackProb = std::exp(trackProb);
    showerProb = std::exp(showerProb);

    return trackProb / (trackProb + showerProb);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NnPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    const float trackLikelihood = GetTrackLikelihood(pCluster);
    const bool isTrackLike(trackLikelihood >= 0.5);

    if (isTrackLike)
        std::cout << "Track-like cluster found: " << trackLikelihood << std::endl;
    else
        std::cout << "Shower-like cluster found: " << (1.0f - trackLikelihood) << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NnPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const /*pPfo*/) const
{
	throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NnPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
