/**
 *  @file   larpandoracontent/LArHelpers/LArTpcGeometryHelper.cc
 *
 *  @brief  Implementation of the geometry helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArTpcGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

LArTpcGeometryHelper *LArTpcGeometryHelper::m_instance = nullptr;

const LArTpcGeometryHelper &LArTpcGeometryHelper::GetInstance(const Algorithm *const pAlgorithm, const LArTransformationPlugin *const pTransform)
{
    if (!m_instance)
        m_instance = new LArTpcGeometryHelper(pAlgorithm, pTransform);

    return *m_instance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TpcHitVolume &LArTpcGeometryHelper::GetTpcHitVolume(const VolumeId &id)
{
    if (m_volumeToTpcMap.find(id) != m_volumeToTpcMap.end())
    {
        return m_volumeToTpcMap.at(id);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArTpcGeometryHelper::LArTpcGeometryHelper(const Algorithm *const pAlgorithm, const LArTransformationPlugin *const pTransform)
{
    const CartesianVector length(359.415,600.019,232.39);
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 0), TpcHitVolume(pAlgorithm, 0, 0, 0, CartesianVector(-182.955,-303.914,115.319),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 2), TpcHitVolume(pAlgorithm, 0, 0, 2, CartesianVector(-182.955,+303.914,115.319),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 4), TpcHitVolume(pAlgorithm, 0, 0, 4, CartesianVector(-182.955,-303.914,347.709),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 6), TpcHitVolume(pAlgorithm, 0, 0, 6, CartesianVector(-182.955,+303.914,347.709),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 8), TpcHitVolume(pAlgorithm, 0, 0, 8, CartesianVector(-182.955,-303.914,580.099),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 10), TpcHitVolume(pAlgorithm, 0, 0, 10, CartesianVector(-182.955,+303.914,580.099),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 12), TpcHitVolume(pAlgorithm, 0, 0, 12, CartesianVector(-182.955,-303.914,812.489),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 14), TpcHitVolume(pAlgorithm, 0, 0, 14, CartesianVector(-182.955,+303.914,812.489),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 16), TpcHitVolume(pAlgorithm, 0, 0, 16, CartesianVector(-182.955,-303.914,1044.88),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 18), TpcHitVolume(pAlgorithm, 0, 0, 18, CartesianVector(-182.955,+303.914,1044.88),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 20), TpcHitVolume(pAlgorithm, 0, 0, 20, CartesianVector(-182.955,-303.914,1277.27),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 0, 22), TpcHitVolume(pAlgorithm, 0, 0, 22, CartesianVector(-182.955,+303.914,1277.27),
        length, pTransform)));

    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 1), TpcHitVolume(pAlgorithm, 0, 1, 1, CartesianVector(+182.955,-303.914,115.319),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 3), TpcHitVolume(pAlgorithm, 0, 1, 3, CartesianVector(+182.955,+303.914,115.319),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 5), TpcHitVolume(pAlgorithm, 0, 1, 5, CartesianVector(+182.955,-303.914,347.709),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 7), TpcHitVolume(pAlgorithm, 0, 1, 7, CartesianVector(+182.955,+303.914,347.709),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 9), TpcHitVolume(pAlgorithm, 0, 1, 9, CartesianVector(+182.955,-303.914,580.099),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 11), TpcHitVolume(pAlgorithm, 0, 1, 11, CartesianVector(+182.955,+303.914,580.099),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 13), TpcHitVolume(pAlgorithm, 0, 1, 13, CartesianVector(+182.955,-303.914,812.489),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 15), TpcHitVolume(pAlgorithm, 0, 1, 15, CartesianVector(+182.955,+303.914,812.489),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 17), TpcHitVolume(pAlgorithm, 0, 1, 17, CartesianVector(+182.955,-303.914,1044.88),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 19), TpcHitVolume(pAlgorithm, 0, 1, 19, CartesianVector(+182.955,+303.914,1044.88),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 21), TpcHitVolume(pAlgorithm, 0, 1, 21, CartesianVector(+182.955,-303.914,1277.27),
        length, pTransform)));
    m_volumeToTpcMap.insert(std::make_pair(VolumeId(0, 1, 23), TpcHitVolume(pAlgorithm, 0, 1, 23, CartesianVector(+182.955,+303.914,1277.27),
         length, pTransform)));
}

} // namespace lar_content

