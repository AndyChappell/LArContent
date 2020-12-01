/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlBranchingAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_BRANCHING_ALGORITHM_H
#define LAR_DL_BRANCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlBranchingAlgorithm class
 */
class DlBranchingAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlBranchingAlgorithm();

    virtual ~DlBranchingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_trackAlgorithms;  ///< The names of the algorithms to run in the track branch
    pandora::StringVector   m_showerAlgorithms; ///< The names of the algortihms to run in the shower branch
    std::string m_trackClusterListName;         ///< The name of the cluster list for track-like clusters
    std::string m_showerClusterListName;        ///< The name of the cluster list for shower-like clusters
    std::string m_outputClusterListName;        ///< The name of the output cluster list
};

} // namespace lar_dl_content

#endif // LAR_DL_BRANCHING_ALGORITHM_H
