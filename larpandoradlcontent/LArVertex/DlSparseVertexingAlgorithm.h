/**
 *  @file   larpandoradlcontent/LArVertexing/DlSparseVertexingAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SPARSE_VERTEXING_ALGORITHM_H
#define LAR_DL_SPARSE_VERTEXING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArObjects/VertexTuple.h"
#include "larpandoradlcontent/LArVertex/DlVertexingBaseAlgorithm.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningTrackShowerIdAlgorithm class
 */
class DlSparseVertexingAlgorithm : public DlVertexingBaseAlgorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DlSparseVertexingAlgorithm();

    ~DlSparseVertexingAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Creates a root ntuple for the training sample to use for determining the vertex location
     *
     *  @return The StatusCode resulting from the function
     */
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief  Estimate the vertex location for an event
     *
     *  @return The StatusCode resulting from the function
     */
    pandora::StatusCode Infer();

    /**
     *  @brief  Determine the region of interest for the hits in a calo hit list
     *
     *  @param  caloHitList The calo hit list to use
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     *  @param  padding The amount of padding to add around the identified region
     */
    void GetHitRegion(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax, const float padding=0.f) const;

    /**
     *  @brief  Get the class from the distance
     *
     *  @param  distance The distance to use
     *
     *  @return The class corresponding to the distance
     */
    int GetClassFromDistance(float distance) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< The name of the calo hit list to use
    std::string m_rootFileName; ///< The name of the root file to write the ntuple to
    std::string m_rootTreeName; ///< The name of the root tree to write the ntuple to
    bool m_trainingMode; ///< Whether the algorithm is in training mode or not
    std::vector<double> m_thresholds; ///< The thresholds for distance classes
};

} // namespace lar_dl_content

#endif // LAR_DL_SPARSE_VERTEXING_ALGORITHM_H
