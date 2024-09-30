/**
 *  @file   larpandoradlcontent/LArVertexing/DlSageVertexingAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SAGE_VERTEXING_ALGORITHM_H
#define LAR_DL_SAGE_VERTEXING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningTrackShowerIdAlgorithm class
 */
class DlSageVertexingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DlSageVertexingAlgorithm();

    ~DlSageVertexingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    bool m_trainingMode;                      ///< Training mode
    std::string m_graphFileName;              ///< Output file name for storing graph information
    std::string m_edgeTreeName;               ///< The name of the tree for storing edge information
    std::string m_nodeTreeName;               ///< The name of the tree for storing node information
    pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    LArDLHelper::TorchModel m_modelU;         ///< The model for the U view
    LArDLHelper::TorchModel m_modelV;         ///< The model for the V view
    LArDLHelper::TorchModel m_modelW;         ///< The model for the W view
    int m_event;                              ///< The current event number
    int m_nClasses;                           ///< The number of distance classes
    std::mt19937 m_rng;                       ///< The random number generator
    std::vector<double> m_thresholds;         ///< Distance class thresholds
};

} // namespace lar_dl_content

#endif // LAR_DL_SAGE_VERTEXING_ALGORITHM_H
