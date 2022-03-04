/**
 *  @file   larpandoracontent/LArUtility/TpcHitVolume.h
 *
 *  @brief  Header file for the list merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TPC_HIT_VOLUME_H
#define LAR_TPC_HIT_VOLUME_H 1

#include "Objects/CaloHit.h"
#include "Objects/CartesianVector.h"
#include "Plugins/LArTransformationPlugin.h"

namespace lar_content
{

/**
 *  @brief  TpcHitVolume class
 */
class TpcHitVolume
{
public:
    /**
     *  @brief Constructor
     *
     *  @param  cryostat The cryostat containing the TPC child volume
     *  @param  tpc The parent TPC volume
     *  @param  child The TPC child volume
     *  @param  center The center of the child TPC volume (x, y, z)
     *  @param  length The full lengths of the sizes of the TPC volume (x, y, z)
     *  @param  pTransform The transformation plugin
     */
    TpcHitVolume(const unsigned int cryostat, const unsigned int tpc, const unsigned int child, const pandora::CartesianVector &center,
        const pandora::CartesianVector &length, const pandora::LArTransformationPlugin *const pTransform);

    TpcHitVolume(const TpcHitVolume &original) = default;

    /**
     *  @brief Destructor
     */
    virtual ~TpcHitVolume();

    /**
     *  @brief  Add (LAr)CaloHits to this child volume. Hits will only be added if they belong to the child volume.
     *
     *  @param  caloHitList The (LAr)CaloHits to be considered for addition to the volume
     */
    void Add(const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Add a (LAr)CaloHit to this child volume. The hit will only be added if it belongs to the child volume.
     *
     *  @param  pCaloHit The (LAr)CaloHit to be considered for addition to the volume
     */
    void Add(const pandora::CaloHit *const pCaloHit);

private:
    /**
     *  @brief  Check if a (LAr)CaloHit is contained within this child volume.
     *
     *  @param  pCaloHit The (LAr)CaloHit to be considered
     *
     *  @return true if the hit is within the volume, false otherwise
     */
    bool Contains(const pandora::CaloHit *const pCaloHit) const;

    /**
     *  @brief  Initialise a view to hits map for a given view.
     *
     *  @param  view The view whose map is to be initialised
     */
    void InitialiseView(const pandora::HitType view);

    typedef std::map<pandora::HitType, pandora::CaloHitList> ViewToHitMap;
    const unsigned int m_cryostat;      ///< The cryostat containing the TPC child volume
    const unsigned int m_tpc;           ///< The parent TPC volume
    const unsigned int m_child;         ///< The TPC child volume
    pandora::CartesianVector m_min;     ///< The minimum coordinate of the TPC child volume XW plane
    pandora::CartesianVector m_max;     ///< The maximum coordinate of the TPC child volume in XW plane
    pandora::CartesianVector m_uMin;     ///< The minimum coordinate of the TPC child volume in XU plane
    pandora::CartesianVector m_uMax;     ///< The maximum coordinate of the TPC child volume in XU plane
    pandora::CartesianVector m_vMin;     ///< The minimum coordinate of the TPC child volume in XV plane
    pandora::CartesianVector m_vMax;     ///< The maximum coordinate of the TPC child volume in XV plane

    ViewToHitMap m_viewToCaloHitMap;    ///< The map from the view to the hits contained within the child volume
};

} // namespace lar_content

#endif // #ifndef LAR_TPC_HIT_VOLUME_H
