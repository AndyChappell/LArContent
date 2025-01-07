/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_KALMAN_CLUSTER_CREATION_ALGORITHM_H
#define LAR_KALMAN_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  KalmanClusterCreationAlgorithm class
 */
class KalmanClusterCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    KalmanClusterCreationAlgorithm();

private:
    typedef std::map<pandora::HitType, pandora::OrderedCaloHitList> ViewOrderedHitsMap;

    /**
     *  @brief  HitAssociation class
     */
    class HitAssociation
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pPrimaryTarget address of the primary target hit
         *  @param  primaryDistanceSquared distance to the primary target hit squared
         */
        HitAssociation(const pandora::CaloHit *const pPrimaryTarget, const float primaryDistanceSquared);

        /**
         *  @brief  Get the primary target
         *
         *  @return the target distance
         */
        const pandora::CaloHit *GetPrimaryTarget() const;

        /**
         *  @brief  Get the primary distance squared
         *
         *  @return the primary distance squared
         */
        float GetPrimaryDistanceSquared() const;

    private:
        const pandora::CaloHit *m_pPrimaryTarget;   ///< the primary target
        float m_primaryDistanceSquared;             ///< the primary distance squared
    };

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Filter out hits below any MIP threshold and organise into an ordered calo hit list
     *
     *  @param  pCaloHitList input hit list
     *  @param  selectedCaloHitList the output list of selected hits
     */
    pandora::StatusCode FilterCaloHits(const pandora::CaloHitList *const pCaloHitList, pandora::OrderedCaloHitList &selectedCaloHitList) const;

    float m_minMipFraction; ///< Minimum fraction of a MIP to consider a hit
    pandora::StringVector m_caloHitListNames; ///< The names of the calo hit lists to cluster
    ViewOrderedHitsMap m_viewHitsMap; ///< Map from the view to the corresponding ordered calo hits
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline KalmanClusterCreationAlgorithm::HitAssociation::HitAssociation(const pandora::CaloHit *const pPrimaryTarget, const float primaryDistanceSquared) :
    m_pPrimaryTarget(pPrimaryTarget),
    m_primaryDistanceSquared(primaryDistanceSquared)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit *KalmanClusterCreationAlgorithm::HitAssociation::GetPrimaryTarget() const
{
    return m_pPrimaryTarget;
}

inline float KalmanClusterCreationAlgorithm::HitAssociation::GetPrimaryDistanceSquared() const
{
    return m_primaryDistanceSquared;
}

} // namespace lar_content

#endif // #ifndef LAR_KALMAN_CLUSTER_CREATION_ALGORITHM_H
