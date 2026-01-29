/**
 *  @file   larpandoracontent/LArMonitoring/ClusterValidationAlgorithm.h
 *
 *  @brief  Header file for the secondary validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_VALIDATION_ALGORITHM_H
#define LAR_CLUSTER_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

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

    virtual ~ClusterValidationAlgorithm();

private:
typedef std::map<const pandora::MCParticle*, int> GenerationMap;
typedef std::map<const pandora::MCParticle*, std::vector<std::pair<const pandora::Cluster*, int>>> ClusterMatchMap;


    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void BuildVisibleHierarchy(const pandora::MCParticle *const pMCParent, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes,
        const int childTier, GenerationMap &generationMap);

    void MatchViewClusters(const GenerationMap &generationMap, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes,
        const pandora::HitType hitType, ClusterMatchMap &clusterMatchMap);

    void FillTree(const GenerationMap &generationMap, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes, 
        const ClusterMatchMap &clusterMatchMapU, const ClusterMatchMap &clusterMatchMapV, const ClusterMatchMap &clusterMatchMapW, 
        const float nuVertexAccuracy, const pandora::MCParticle *const pMCNu);

    void GetHitsOfType(const pandora::CaloHitList &inputList, const pandora::HitType hitType, pandora::CaloHitList &outputList);

    void GetMatchingMetrics(const pandora::CaloHitList &targetHitList, const std::vector<std::pair<const pandora::Cluster*, int>> matches,
         pandora::FloatVector &completeness, pandora::FloatVector &purity, pandora::IntVector &nSharedHits, pandora::IntVector &isReconstructableView);

    bool FindMCNode(const pandora::MCParticle *const pMCParticle, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyMCNodes, 
                    const LArHierarchyHelper::MCHierarchy::Node *&pMCNode);

    int m_eventNumber;
    bool m_writeFile;
    std::string m_fileName;
    std::string m_treeName;
    unsigned int m_minRecoHits;        ///< Minimum number of reconstructed primary good hits
    unsigned int m_minRecoHitsPerView; ///< Minimum number of reconstructed hits for a good view
    unsigned int m_minRecoGoodViews;   ///< Minimum number of reconstructed primary good views
    bool m_removeRecoNeutrons;         ///< Whether to remove reconstructed neutrons and their downstream particles
};

} // namespace lar_content

#endif // LAR_CLUSTER_VALIDATION_ALGORITHM_H
