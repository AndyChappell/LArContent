/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include <cmath>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoradlcontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlPfoCharacterisationAlgorithm::DlPfoCharacterisationAlgorithm() 
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
            const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
            const float pTrack{pLArCaloHit->GetTrackProbability()};
            const float pShower{pLArCaloHit->GetShowerProbability()};
            trackProb += pTrack;
            showerProb += pShower;
        }
    }

    return trackProb / (trackProb + showerProb);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    const float trackLikelihood = GetTrackLikelihood(pCluster);
    const bool isTrackLike(trackLikelihood >= 0.5);

    return isTrackLike;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const /*pPfo*/) const
{
	throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content

