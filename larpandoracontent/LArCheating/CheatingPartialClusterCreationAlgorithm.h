/**
 *  @file   larpandoracontent/LArCheating/CheatingPartialClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_PARTIAL_CLUSTER_CREATION_ALGORITHM_H
#define LAR_CHEATING_PARTIAL_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingPartialClusterCreationAlgorithm class
 */
class CheatingPartialClusterCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingPartialClusterCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::MCParticle *, pandora::CaloHitList> MCParticleToHitListMap;

    /**
     *  @brief  Create map between each (primary) MC particle and associated calo hits
     *
     *  @param  mcParticleToHitListMap to receive the mc particle to hit list map
     */
    void GetMCParticleToHitListMap(MCParticleToHitListMap &mcParticleToHitListMap) const;

    /**
     *  @brief  Create clusters based on information in the mc particle to hit list map
     *
     *  @param  mcParticleToHitListMap the mc particle to hit list map
     */
    void CreateClusters(const MCParticleToHitListMap &mcParticleToHitListMap) const;

    pandora::IntVector m_pdgExclusionVector;    ///< The PDG codes to avoid clustering
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_PARTIAL_CLUSTER_CREATION_ALGORITHM_H
