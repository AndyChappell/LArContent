/**
 *  @file   larpandoracontent/LArHelpers/LArCalorimetryHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArCalorimetryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

void LArCalorimetryHelper::GetPrincipalAxes(const CaloHitList &caloHitList, CartesianVector &origin, CartesianVector &longitudinal,
    CartesianVector &transverse)

{
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(caloHitList, origin, eigenValues, eigenVectors);

    longitudinal = eigenVectors[0];
    transverse = eigenVectors[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector LArCalorimetryHelper::GetLongitudinalAdcProfile(const CaloHitList &caloHitList, const CartesianVector &origin, const CartesianVector &axis,
    const float binWidth, const float interactionLength)
{
    FloatVector projVector, adcVector;
    float min{std::numeric_limits<float>::max()}, max{-std::numeric_limits<float>::max()};
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float proj{axis.GetDotProduct(pCaloHit->GetPositionVector() - origin)};
        min = std::min(min, proj);
        max = std::max(max, proj);
        projVector.emplace_back(proj);
        adcVector.emplace_back(pCaloHit->GetInputEnergy());
    }
    if (adcVector.empty())
        return FloatVector();
    FloatVector sortedAdcVector(adcVector);
    std::sort(sortedAdcVector.begin(), sortedAdcVector.end());
    float medianAdc{1.f};
    if (sortedAdcVector.size() % 2)
    {
        const size_t midPoint{sortedAdcVector.size() / 2};
        if (sortedAdcVector[midPoint] > std::numeric_limits<float>::epsilon())
            medianAdc = sortedAdcVector[midPoint];
    }
    else
    {
        const size_t a{sortedAdcVector.size() / 2 - 1};
        const size_t b{sortedAdcVector.size() / 2};
        if (sortedAdcVector[a] > std::numeric_limits<float>::epsilon() && sortedAdcVector[b] > std::numeric_limits<float>::epsilon())
            medianAdc = 0.5f * (sortedAdcVector[a] + sortedAdcVector[b]);
    }

    // Adjust min and max to avoid underflow/overflow issues
    min -= std::numeric_limits<float>::epsilon();
    max += std::numeric_limits<float>::epsilon();
    float binSize{binWidth * interactionLength};
    const size_t N{static_cast<size_t>(std::ceil((max - min) / binSize))};
    if (N < 1)
        return FloatVector();

    // Make ADC profile
    FloatVector profile(N);
    for (size_t i = 0; i < projVector.size(); ++i)
    {
        const int bin{static_cast<int>((projVector[i] - min) / binSize)};
        profile[bin] += adcVector[i] / medianAdc;
    }

    return profile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector LArCalorimetryHelper::GetTransverseAdcProfile(const CaloHitList &caloHitList, const CartesianVector &origin, const CartesianVector &axis,
    const float binWidth, const float interactionLength)
{
    FloatVector projVector, adcVector;
    float min{std::numeric_limits<float>::max()}, max{-std::numeric_limits<float>::max()};
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float proj{axis.GetDotProduct(pCaloHit->GetPositionVector() - origin)};
        min = std::min(min, proj);
        max = std::max(max, proj);
        projVector.emplace_back(proj);
        adcVector.emplace_back(pCaloHit->GetInputEnergy());
    }
    // Make the bounds symmetric
    if (std::abs(min) > std::abs(max))
        max = -min;
    else
        min = -max;
    if (adcVector.empty())
        return FloatVector();
    FloatVector sortedAdcVector(adcVector);
    std::sort(sortedAdcVector.begin(), sortedAdcVector.end());
    float medianAdc{1.f};
    if (sortedAdcVector.size() % 2)
    {
        const size_t midPoint{sortedAdcVector.size() / 2};
        if (sortedAdcVector[midPoint] > std::numeric_limits<float>::epsilon())
            medianAdc = sortedAdcVector[midPoint];
    }
    else
    {
        const size_t a{sortedAdcVector.size() / 2 - 1};
        const size_t b{sortedAdcVector.size() / 2};
        if (sortedAdcVector[a] > std::numeric_limits<float>::epsilon() && sortedAdcVector[b] > std::numeric_limits<float>::epsilon())
            medianAdc = 0.5f * (sortedAdcVector[a] + sortedAdcVector[b]);
    }

    // Adjust min and max to avoid underflow/overflow issues
    min -= std::numeric_limits<float>::epsilon();
    max += std::numeric_limits<float>::epsilon();
    float binSize{binWidth * interactionLength};
    size_t N{static_cast<size_t>(std::ceil((max - min) / binSize))};
    if (N < 1)
    {
        return FloatVector();
    }
    else if ((N % 2) == 0)
    {
        // Force symmetric binning
        min -= 0.5f * binSize;
        max += 0.5f * binSize;
        N = static_cast<size_t>(std::ceil((max - min) / binSize));
    }

    // Make ADC profile
    FloatVector profile(N);
    for (size_t i = 0; i < projVector.size(); ++i)
    {
        const int bin{static_cast<int>((projVector[i] - min) / binSize)};
        profile[bin] += adcVector[i] / medianAdc;
    }

    return profile;
}

} // namespace lar_content
