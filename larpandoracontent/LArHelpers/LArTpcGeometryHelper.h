/**
 *  @file   larpandoracontent/LArHelpers/LArTpcGeometryHelper.h
 *
 *  @brief  Header file for the geometry helper class.
 *
 *  $Log: $
 */
#ifndef LAR_TPC_GEOMETRY_HELPER_H
#define LAR_TPC_GEOMETRY_HELPER_H 1

#include "Plugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArUtility/TpcHitVolume.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  LArTpcGeometryHelper class
 */
class LArTpcGeometryHelper
{
public:
    typedef std::tuple<unsigned int, unsigned int, unsigned int> VolumeId;

    /**
     *  @brief  Get the instance of the TPC geometry helper.
     *
     *  @param  pTransform The transformation plugin
     *
     *  @return The TPC geometry helper
     */
    static const LArTpcGeometryHelper &GetInstance(const pandora::LArTransformationPlugin *const pTransform);

    /**
     *  @brief  Retrieve a TpcHitVolume.
     *
     *  @param  id The volume id for the TPC to be retrieved
     *
     *  @return The TPC requested
     */
    TpcHitVolume &GetTpcHitVolume(const VolumeId &id);

private:
    typedef std::map<VolumeId, TpcHitVolume> VolumeIdToTpcMap;

    /**
     *  Constructor
     *
     *  @param  pTransform The transformation plugin
     */
    LArTpcGeometryHelper(const pandora::LArTransformationPlugin *const pTransform);

    static LArTpcGeometryHelper *m_instance;    ///< Singleton describing the TPC geometry
    VolumeIdToTpcMap m_volumeToTpcMap;          ///< The map from ID to TPC
};

} // namespace lar_content

#endif // #ifndef LAR_TPC_GEOMETRY_HELPER_H

