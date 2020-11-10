/**
 *  @file   larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.h
 *
 *  @brief  Header file for the cluster validation class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_VALIDATION_ALGORITHM_H
#define LAR_CLUSTER_VALIDATION_ALGORITHM_H 1

#include "Objects/CaloHit.h"
#include "Objects/MCParticle.h"

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

#include <functional>

namespace lar_content
{

/**
 *  @brief  ClusterValidationAlgorithm class
 */
class ClusterValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ClusterValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~ClusterValidationAlgorithm();

protected:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    /**
     *  @brief  Determine the mapping between MC particles and calo hits
     *
     *  @param  pMCParticleList The input MC particle list
     *  @param  pCaloHitList The input calo hit list
     *  @param  targetMCToHitsMap The ouput map between selected MC particles and calo hits
     *  @param  allMCToHitsMap The output map between all MC particles and calo hits
     */
    void GetMCToHitsMaps(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList,
        LArMCParticleHelper::MCContributionMap &targetMCToHitsMap, LArMCParticleHelper::MCContributionMap &allMCToHitsMap);

    /**
     *  @brief  Process the MC to calo hit mapping and output the result
     *
     *  @param  targetMCToHitsMap The input map between selected MC particles and calo hits
     *  @param  clusterToHitsMap The input map between clusters and calo hits
     */
    void ProcessOutput(const LArMCParticleHelper::MCContributionMap &targetMCToHitsMap,
        const LArMCParticleHelper::ClusterContributionMap &clusterToHitsMap);

    typedef std::unordered_map<const pandora::Cluster*, unsigned int> ClusterToIdMap;

    std::string                                             m_caloHitListName;              ///< Name of input calo hit list
    std::string                                             m_mcParticleListName;           ///< Name of input MC particle list
    pandora::StringVector                                   m_clusterListNames;             ///< Names of input cluster lists

    bool                                                    m_writeToTree;                  ///< Whether to write validation details to tree
    std::string                                             m_treeName;                     ///< Name of output tree
    std::string                                             m_fileName;                     ///< Name of output file

    int                                                     m_eventNumber;                  ///< The event number

    std::function<bool(const pandora::MCParticle *const)>   m_criteria;                     ///< The selection criteria for MC particles
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_HELPER_H

