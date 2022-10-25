/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/HitCorrelationAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_CORRELATION_ALGORITHM_H
#define LAR_HIT_CORRELATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  HitCorrelationAlgorithm class
 */
class HitCorrelationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HitCorrelationAlgorithm();

    ~HitCorrelationAlgorithm();

private:
    class TpcChildVolume;
    typedef std::map<const pandora::CaloHit *, pandora::CaloHitList> HitMap;

    /**
     *  @brief  Find pairs of 2D hits that might originate from the same 3D hit
     *
     *  @param  caloHitList1 A list of CaloHits in a view, sorted by x position
     *  @param  caloHitList2 A second  list of CaloHits in a view, sorted by x position
     *  @param  hitMap A map from CaloHits in one view to CaloHits in the other
     **/
    void Correlate(const pandora::CaloHitList &caloHitList1, const pandora::CaloHitList &caloHitList2, HitMap &hitMap) const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< The name of the calo hit list to use (this should be a multi-HitType list)
    std::map<unsigned int, TpcChildVolume> m_volumeMap; ///< A map containing the various TPC child volumes
    std::map<const pandora::CaloHit *, bool> m_availabilityMap; ///< A map indicating the availability of hits

    class TpcChildVolume
    {
    public:
        /**
         *  @brief  Default constructor
         */
        TpcChildVolume() = default;

        /**
         *  @brief  Add a CaloHit to the TPC child volume
         */
        void AddCaloHit(const pandora::CaloHit *pCaloHit);

        /**
         *  @brief  Retrieve the collection of hits of a given type in this TPC child volume
         *
         *  @param  type The HitType of the hits to retrieve
         *
         *  @return The list of hits
         */
        const pandora::CaloHitList GetCaloHits(const pandora::HitType type) const;

    private:
        pandora::CaloHitList m_hitsU;
        pandora::CaloHitList m_hitsV;
        pandora::CaloHitList m_hitsW;
    };
};

} // namespace lar_content

#endif // #ifndef LAR_HIT_CORRELATION_ALGORITHM_H
