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

    /**
     *  @brief  Comparator for CartesianVector to be used in std::map
     */
    struct VertexComparator
    {
        bool operator()(const pandora::CartesianVector &a, const pandora::CartesianVector &b) const
        {
            if (a.GetX() != b.GetX())
                return a.GetX() < b.GetX();
            else
                return a.GetZ() < b.GetZ();
        }
    };
    typedef std::map<const pandora::CartesianVector, pandora::CaloHitList, VertexComparator> VertexHitsMap;
    typedef std::map<const pandora::CaloHit*, int> CondensationPointMap;
    typedef std::map<const pandora::CaloHit*, pandora::IntVector> HitVertexLabelMap;
    typedef std::map<const pandora::CaloHit*, pandora::FloatVector> HitVertexWeightMap;

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

    /**
     *  @brief  Consolidate any duplicate vertices and associate the relevant hits to the consolidated vertex
     *
     *  @param  mcToMatchedVertexMap the mapping between MC particles and their matched vertex positions
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     *  @param  vertexHitsMap the output mapping between consolidated vertex positions and their associated hits
     */
    void ConsolidateVertices(const MCVertexMap &mcToMatchedVertexMap, const LArMCParticleHelper::MCContributionMap &mcToHitsMap, VertexHitsMap &vertexHitsMap) const;

    /**
     *  @brief  Filter out particles that are neutron-induced
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     */
    void FilterNeutronInducedParticles(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    /**
     *  @brief  Make the training file from the hits and consolidated vertices
     *
     *  @param  treeName the name of the ROOT tree
     *  @param  vertexHitsMap the mapping between consolidated vertex positions and their associated hits
     *  @param  mcToMatchedVertexMap the mapping between MC particles and their matched vertex positions
     *  @param  fullCaloHitList the full list of calo hits from which the training sample is drawn
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     */
    void MakeTrainingFile(const std::string &treeName, const VertexHitsMap &vertexHitsMap, const MCVertexMap &mcToMatchedVertexMap,
        const pandora::CaloHitList &fullCaloHitList, const LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    /**
     *  @brief  Visualize the hits and vertices for each MC particle
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     *  @param  mcVertexMap the mapping between MC particles and their projected true vertex positions
     *  @param  mcToMatchedVertexMap the output mapping between MC particles and their matched vertex positions
     */
    void VisualizeByMC(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCVertexMap &mcVertexMap, const MCVertexMap &mcToMatchedVertexMap) const;

    /**
     *  @brief  Visualize the consolidated vertices and their associated hits
     *
     *  @param  vertexHitsMap the mapping between consolidated vertex positions and their associated hits
     */
    void VisualizeByVertex(const VertexHitsMap &vertexHitsMap) const;

    bool m_trainingMode;   ///< Whether or not the algorithm is in training mode
    bool m_visualize;           ///< Whether or not to visualize the candidate vertices
    pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    std::string m_rootTreeName; ///< The ROOT tree name
    std::string m_rootFileName; ///< The ROOT file name
    std::mt19937 m_rng;         ///< The random number generator
};

} // namespace lar_dl_content

#endif // LAR_DL_VERTEX_CONDENSATION_ALGORITHM_H
