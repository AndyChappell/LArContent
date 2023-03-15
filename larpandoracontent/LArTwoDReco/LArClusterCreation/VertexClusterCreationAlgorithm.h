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

#include <unordered_map>

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
    void IdentifyTrackStubs(const pandora::CaloHitList &caloHitList, const pandora::Vertex &vertex) const;
    bool ClusterHits(const pandora::CaloHit *const pSeedHit, const pandora::CaloHitVector &caloHitVector, const pandora::CartesianVector &vertex) const;
    float GetClosestApproach(const pandora::CaloHitVector &caloHitVector, const pandora::CartesianVector &vertex) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< The name of the calo hit list
    std::string m_vertexListName;   ///< The name of the vertex list
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_CLUSTER_CREATION_ALGORITHM_H
