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
#include "Pandora/Algorithm.h"
#include "Plugins/LArTransformationPlugin.h"

namespace lar_content
{

/**
 *  @brief  TpcHitVolume class
 */
class TpcHitVolume
{
public:
    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterMap;

    /**
     *  @brief Constructor
     *
     *  @param  pAlgorithm The algorithm creating the object
     *  @param  cryostat The cryostat containing the TPC child volume
     *  @param  tpc The parent TPC volume
     *  @param  child The TPC child volume
     *  @param  center The center of the child TPC volume (x, y, z)
     *  @param  length The full lengths of the sizes of the TPC volume (x, y, z)
     *  @param  pTransform The transformation plugin
     */
    TpcHitVolume(const pandora::Algorithm *const pAlgorithm, const unsigned int cryostat, const unsigned int tpc, const unsigned int child,
        const pandora::CartesianVector &center, const pandora::CartesianVector &length, const pandora::LArTransformationPlugin *const pTransform);

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

    /**
     *  @brief  Get the relationships between clusters in the different views.
     *
     *  @param  signalMap The output map from clusters in one view to clusters in another that constitute signal
     *  @param  backgroundMap The output map from clusters in one view to clusters in another that constiture background
     */
    void GetCombinatorics(ClusterMap &signalMap, ClusterMap &backgroundMap) const;

    /**
     *  @brief  Retrieve the hits from a given view in this volume
     *
     *  @param  view The view from which hits should be retrieved
     *  @param  caloHitList The output CaloHitList
     */
    void GetHitList(const pandora::HitType view, pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Take all of the hits from a given view in this volume and transforms them to the local coorindate system
     *
     *  @param  view The view for which coordinates should be produced
     *  @param  localCoords The output vector of coordinates
     */
    void GetLocalCoordinates(const pandora::HitType view, pandora::CartesianPointVector &localCoords) const;

    /**
     *  @brief  Take a given cluster and transform its hits into the local coorindate system
     *
     *  @param  view The view for which coordinates should be produced
     *  @param  localCoords The output vector of coordinates
     */
    void GetLocalCoordinates(const pandora::Cluster *const pCluster, pandora::CartesianPointVector &localCoords) const;

    /**
     *  @brief  Take a given cluster and transform its hits into the local coorindate system
     *
     *  @param  view The view for which coordinates should be produced
     *  @param  xMin The minimum x coordinate for hits to be transformed
     *  @param  xMax The maximum x coordinate for hits to be transformed
     *  @param  localCoords The output vector of coordinates
     */
    void GetLocalCoordinates(const pandora::Cluster *const pCluster, const float xMin, const float xMax,
        pandora::CartesianPointVector &localCoords) const;

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
     *  @brief  Transform the hit coordinate to the local coorindate system
     *
     *  @param  pCaloHit The (LAr)CaloHit to transform
     *  @param  xMin The minimum x coordinate for hits to be transformed
     *  @param  xMax The maximum x coordinate for hits to be transformed
     *  @param  localCoords The output vector of coordinates in which the transformed coordinate should be stored
     */
    void GetLocalCoordinate(const pandora::CaloHit *pCaloHit, const float xMin, const float xMax, pandora::CartesianPointVector &localCoords) const;

    /**
     *  @brief  Retrieve the main MC particle associated with a cluster, along with the weight this particle contributes to the cluster
     *
     *  @param  pCluster The cluster for which the MC particle should be determined
     *  @param  weight The output weight that the returned particle contributes to the cluster
     *
     *  @return The MCParticle that best represents the cluster
     */
    const pandora::MCParticle *GetMainMCParticle(const pandora::Cluster *const pCluster, float &weight) const;

    /**
     *  @brief  Initialise a view to hits map for a given view.
     *
     *  @param  view The view whose map is to be initialised
     */
    void InitialiseView(const pandora::HitType view);

    /**
     *  @brief  Debug function that prints out various hit and volume information.
     *
     *  @param  pCaloHit The (LAr)CaloHit about which information should be printed
     */
    void _Echo(const pandora::CaloHit *const pCaloHit) const;

    typedef std::map<pandora::HitType, pandora::CaloHitList> ViewToHitMap;

    const pandora::Algorithm *const m_pAlgorithm;   ///< The parent algorithm
    const unsigned int m_cryostat;                  ///< The cryostat containing the TPC child volume
    const unsigned int m_tpc;                       ///< The parent TPC volume
    const unsigned int m_child;                     ///< The TPC child volume
    pandora::CartesianVector m_min;                 ///< The minimum coordinate of the TPC child volume XW plane
    pandora::CartesianVector m_max;                 ///< The maximum coordinate of the TPC child volume in XW plane
    pandora::CartesianVector m_uMin;                ///< The minimum coordinate of the TPC child volume in XU plane
    pandora::CartesianVector m_uMax;                ///< The maximum coordinate of the TPC child volume in XU plane
    pandora::CartesianVector m_vMin;                ///< The minimum coordinate of the TPC child volume in XV plane
    pandora::CartesianVector m_vMax;                ///< The maximum coordinate of the TPC child volume in XV plane
    ViewToHitMap m_viewToCaloHitMap;                ///< The map from the view to the hits contained within the child volume
};

} // namespace lar_content

#endif // #ifndef LAR_TPC_HIT_VOLUME_H
