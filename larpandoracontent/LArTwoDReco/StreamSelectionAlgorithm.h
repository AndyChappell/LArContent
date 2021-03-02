/**
 *  @file   larpandoradlcontent/LArTwoDReco/StreamSelectionAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_STREAM_SELECTION_ALGORITHM_H
#define LAR_STREAM_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  StreamSelectionAlgorithm class
 */
class StreamSelectionAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    StreamSelectionAlgorithm();

    virtual ~StreamSelectionAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_trackClusterListName;     ///< The name of the output cluster list for track-like clusters
    std::string m_showerClusterListName;    ///< The name of the output cluster list for shower-like clusters
    bool        m_useTracksAsOutputList;    ///< Whether or not to use tracks for the output cluster list
};

} // namespace lar_content

#endif // LAR_STREAM_SELECTION_ALGORITHM_H
