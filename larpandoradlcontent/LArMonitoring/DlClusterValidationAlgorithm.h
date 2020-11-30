/**
 *  @file   larpandoradlcontent/LArMonitoring/DlClusterValidationAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id cluster validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_CLUSTER_VALIDATION_ALGORITHM_H
#define LAR_DL_CLUSTER_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlClusterValidationAlgorithm class
 */
class DlClusterValidationAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlClusterValidationAlgorithm();

    virtual ~DlClusterValidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< Name of input calo hit list
    std::string m_outputFileName;   ///< Name of output validation root file
    std::string m_outputTreeName;   ///< Name of output ROOT tree
};

} // namespace lar_dl_content

#endif // LAR_DL_CLUSTER_VALIDATION_ALGORITHM_H
