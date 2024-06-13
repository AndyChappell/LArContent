/**
 *  @file   larpandoracontent/LArTwoDReco/LArAssociatedCluster/DlVertexAssociatedClusterAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_VERTEX_ASSOCIATED_CLUSTER_ALGORITHM_H
#define LAR_DL_VERTEX_ASSOCIATED_CLUSTER_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

using namespace lar_content;

namespace lar_dl_content
{

/**
 *  @brief  DlVertexAssociatedClusterAlgorithm class
 */
class DlVertexAssociatedClusterAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlVertexAssociatedClusterAlgorithm();

    ~DlVertexAssociatedClusterAlgorithm();

private:
    pandora::StatusCode Run();
    void IdentifyAssociatedClusters(const pandora::ClusterList &clusterList, const pandora::VertexList &vertexList) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_clusterListName;      ///< The name of the cluster list
    std::string m_vertexListName;       ///< The name of the vertex list
    float m_rInfluenceSquared;          ///< The square of the radius of influence of a vertex
    float m_halfBoundingEdge;           ///< Half the length of the bounding box edge
    int m_nBins;                        ///< The number of bins along a bounding edge
    bool m_trainingMode;                ///< Whether or not we're running in training mode
    std::string m_rootTreeName;         ///< The ROOT tree name
    std::string m_rootFileName;         ///< The ROOT file name
};

} // namespace lar_content

#endif // #ifndef LAR_DL_VERTEX_ASSOCIATED_CLUSTER_ALGORITHM_H
