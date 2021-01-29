/**
 *  @file   larpandoracontent/LArMonitoring/SimpleClusterValidationAlgorithm.h
 *
 *  @brief  Header file for the simple cluster validation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_CLUSTER_VALIDATION_ALGORITHM_H
#define LAR_SIMPLE_CLUSTER_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  SimpleClusterValidationAlgorithm class
 */
class SimpleClusterValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SimpleClusterValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~SimpleClusterValidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Generates the metrics for cluster validation
     *
     *  @param  pClusterList the input list of clusters
     *  @param  hitToMCMap the input map from calo hits to MC particles
     *  @param  mcToHitMap the input map from MC particles to calo hits
     *
     *  @return The StatusCode
     */
    pandora::StatusCode GenerateMetrics(const pandora::ClusterList *pClusterList, const LArMCParticleHelper::CaloHitToMCMap &hitToMCMap,
        const LArMCParticleHelper::MCContributionMap &mcToHitMap);

    /**
     *  @brief  Get the maps between MC particles and calo hits
     *
     *  @param  pMCList the input list of MC particles
     *  @param  pClusterList the input list of clusters
     *  @param  hitToMCMap the output map from calo hits to MC particles
     *  @param  mcToHitMap the output map from MC particles to calo hits
     *
     *  @return The StatusCode
     */
    pandora::StatusCode GetRecoTruthMaps(const pandora::MCParticleList *pMCList, const pandora::ClusterList *pClusterList,
        LArMCParticleHelper::CaloHitToMCMap &hitToMCMap, LArMCParticleHelper::MCContributionMap &mcToHitMap);

    pandora::StringVector   m_clusterListNames;         ///< The list of cluster list names
    std::string             m_filename;                 ///< The output ROOT filename
    std::string             m_treeName;                 ///< The output ROOT tree name
    bool                    m_writeFile;                ///< Whether or not to write the ROOT tree to file
};

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_CLUSTER_VALIDATION_ALGORITHM_H
