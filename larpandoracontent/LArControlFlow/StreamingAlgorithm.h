/**
 *  @file   larpandoradlcontent/LArTwoDReco/StreamingAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_STREAMING_ALGORITHM_H
#define LAR_STREAMING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  StreamingAlgorithm class
 */
class StreamingAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    StreamingAlgorithm();

    virtual ~StreamingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_trackAlgorithms;  ///< The names of the algorithms to run in the track branch
    pandora::StringVector   m_showerAlgorithms; ///< The names of the algortihms to run in the shower branch
    std::string m_trackClusterListName;         ///< The name of the cluster list for track-like clusters
    std::string m_showerClusterListName;        ///< The name of the cluster list for shower-like clusters
    std::string m_outputClusterListName;        ///< The name of the output cluster list
};

} // namespace lar_content

#endif // LAR_STREAMING_ALGORITHM_H
