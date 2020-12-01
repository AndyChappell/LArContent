/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlClusterStreamingAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_CLUSTER_STREAMING_ALGORITHM_H
#define LAR_DL_CLUSTER_STREAMING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlClusterStreamingAlgorithm class
 */
class DlClusterStreamingAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlClusterStreamingAlgorithm();

    virtual ~DlClusterStreamingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_trackClusterListName;     ///< The name of the output cluster list for track-like clusters
    std::string m_showerClusterListName;    ///< The name of the output cluster list for shower-like clusters
    bool        m_useTracksAsOutputList;    ///< Whether or not to use tracks for the output cluster list
};

} // namespace lar_dl_content

#endif // LAR_DL_CLUSTER_STREAMING_ALGORITHM_H
