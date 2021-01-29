/**
 *  @file   larpandoradlcontent/LArCheating/CheatingClusterStreamingAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_CLUSTER_STREAMING_ALGORITHM_H
#define LAR_CHEATING_CLUSTER_STREAMING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingClusterStreamingAlgorithm class
 */
class CheatingClusterStreamingAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingClusterStreamingAlgorithm();

    virtual ~CheatingClusterStreamingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_trackClusterListName;     ///< The name of the output cluster list for track-like clusters
    std::string m_showerClusterListName;    ///< The name of the output cluster list for shower-like clusters
};

} // namespace lar_content

#endif // LAR_CHEATING_CLUSTER_STREAMING_ALGORITHM_H
