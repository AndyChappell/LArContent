/**
 *  @file   larpandoradlcontent/LArVertexCondensation/DlVertexCondensationAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_VERTEX_CONDENSATION_ALGORITHM_H
#define LAR_DL_VERTEX_CONDENSATION_ALGORITHM_H 1

#include <Eigen/Dense>

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
class DlVertexCondensationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DlVertexCondensationAlgorithm();

    ~DlVertexCondensationAlgorithm();

private:
    typedef std::map<const pandora::MCParticle*, pandora::CartesianVector> MCVertexMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /**
     *  @brief  Get the projected true vertex positions in a given view
     *
     *  @param  pMC the MC particle for which the vertex should be found
     *  @param  view the TPC view
     *  @param  mcVertexMap the output mapping between MC particles and projected true vertex positions
     */
    void GetProjectedTrueVertices(const pandora::MCParticle *const pMC, const pandora::HitType view, MCVertexMap &mcVertexMap) const;

    /**
     *  @brief  Match each vertex to the closest hit and return the distance
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     *  @param  mcVertexMap the mapping between MC particles and their projected true vertex positions
     *  @param  mcToMatchedVertexMap the output mapping between MC particles and their matched vertex positions
     */
    void MatchHitToVertex(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCVertexMap &mcVertexMap, MCVertexMap &mcToMatchedVertexMap) const;

    bool m_trainingMode;   ///< Whether or not the algorithm is in training mode
    bool m_visualize;           ///< Whether or not to visualize the candidate vertices
    pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    std::string m_rootTreeName; ///< The ROOT tree name
    std::string m_rootFileName; ///< The ROOT file name
    std::mt19937 m_rng;         ///< The random number generator
};

} // namespace lar_dl_content

#endif // LAR_DL_VERTEX_CONDENSATION_ALGORITHM_H
