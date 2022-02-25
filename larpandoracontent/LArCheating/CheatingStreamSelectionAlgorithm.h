/**
 *  @file   larpandoracontent/LArTwoDReco/CheatingStreamSelectionAlgorithm.h
 *
 *  @brief  Header file for cheating stream selection algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_STREAM_SELECTION_ALGORITHM_H
#define LAR_CHEATING_STREAM_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/StreamSelectionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingStreamSelectionAlgorithm class
 */
class CheatingStreamSelectionAlgorithm : public lar_content::StreamSelectionAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingStreamSelectionAlgorithm() = default;

    virtual ~CheatingStreamSelectionAlgorithm() = default;

protected:
    /**
     *  @brief  Allocate a cluster to the appropriate streams.
     *
     *  @param  pCluster The cluster to allocate to a stream
     *
     *  @return The StatusCode
     */
    virtual pandora::StatusCode AllocateToStreams(const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_trackListName;  ///< The name of the track list
    std::string m_showerListName; ///< The name of the shower list
};

} // namespace lar_content

#endif // LAR_CHEATING_STREAM_SELECTION_ALGORITHM_H
