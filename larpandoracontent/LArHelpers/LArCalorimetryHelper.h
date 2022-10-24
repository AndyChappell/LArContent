/**
 *  @file   larpandoracontent/LArHelpers/LArCalorimetryHelper.h
 *
 *  @brief  Header file for the cluster helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CALORIMETRY_HELPER_H
#define LAR_CALORIMETRY_HELPER_H 1

#include "Objects/CaloHit.h"

namespace lar_content
{

/**
 *  @brief  LArCalorimetryHelper class
 */
class LArCalorimetryHelper
{
public:
    /**
     *  @brief  Get the principlal axes for a collection of 2D calo hits
     *
     *  @param  caloHitList The collection of hits for which axes should be determined
     *  @param  origin The output origin of the principal axes
     *  @param  longitudinal The output longitudinal axis for the hits
     *  @param  transverse The output transverse axis for the hits
     */
    static void GetPrincipalAxes(const pandora::CaloHitList &caloHitList, pandora::CartesianVector &origin, pandora::CartesianVector &longitudinal,
        pandora::CartesianVector &transverse);

    /**
     *  @brief  Get the energy profile for a collection of 2D calo hits along a given axis
     *
     *  @param  caloHitList The collection of hits for which axes should be determined
     *  @param  origin The origin of the axis along which the profile should be calculated
     *  @param  axis The axis along which the profile should be calculated
     *  @param  binWidth The size of each bin in interaction lengths
     *  @param  interactionLength The interaction length
     *
     *  @return The energy profile binned in interaction lengths
     */
    static pandora::FloatVector GetLongitudinalAdcProfile(const pandora::CaloHitList &caloHitLits, const pandora::CartesianVector &origin,
        const pandora::CartesianVector &axis, const float binWidth, const float interactionLength);

    /**
     *  @brief  Get the transverse energy profile for a collection of 2D calo hits along a given axis
     *
     *  @param  caloHitList The collection of hits for which axes should be determined
     *  @param  origin The origin of the axis along which the profile should be calculated
     *  @param  axis The axis along which the profile should be calculated
     *  @param  binWidth The size of each bin in interaction lengths
     *  @param  interactionLength The interaction length
     *
     *  @return The energy profile binned in interaction lengths, zero padded so the origin is in the central bin
     */
    static pandora::FloatVector GetTransverseAdcProfile(const pandora::CaloHitList &caloHitLits, const pandora::CartesianVector &origin,
        const pandora::CartesianVector &axis, const float binWidth, const float interactionLength);
};

} // namespace lar_content

#endif // #ifndef LAR_CALORIMETRY_HELPER_H
