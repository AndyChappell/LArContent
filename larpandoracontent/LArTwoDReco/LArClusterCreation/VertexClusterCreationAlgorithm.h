/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/VertexClusterCreationAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_CLUSTER_CREATION_ALGORITHM_H
#define LAR_VERTEX_CLUSTER_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include <vector>

namespace lar_content
{

/**
 *  @brief  VertexClusterCreationAlgorithm class
 */
class VertexClusterCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexClusterCreationAlgorithm();

private:
    pandora::StatusCode Run();
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void IdentifyTrackStubs(const pandora::CaloHitList &caloHitList, const pandora::Vertex &vertex) const;
    void IdentifyConnectingPaths(const pandora::CaloHitList &caloHitList, const pandora::Vertex &vertex) const;
    void ClusterHits(const pandora::CaloHitVector &caloHitVector) const;
    float GetClosestApproach(const pandora::CaloHitVector &caloHitVector, const pandora::CartesianVector &vertex) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< The name of the calo hit list
    std::string m_vertexListName;   ///< The name of the vertex list

    float m_hitRadii; ///< The search radius to collect hits around the vertex
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_CLUSTER_CREATION_ALGORITHM_H
