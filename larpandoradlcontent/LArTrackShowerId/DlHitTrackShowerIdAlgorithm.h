/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

namespace lar_dl_content
{

/**
 *  @brief  DlHitTrackShowerIdAlgorithm class
 */
class DlHitTrackShowerIdAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlHitTrackShowerIdAlgorithm();

    virtual ~DlHitTrackShowerIdAlgorithm();

private:
    typedef std::unordered_map<const pandora::MCParticle *, const pandora::MCParticle *> MCFoldingMap;

    enum Category
    {
        UNINITIALISED = 0,
        MIP = 1,
        HIP = 2,
        SHOWER = 3,
        LOW_E = 4, // Michels and Deltas
        DIFFUSE = 5 // Compton scattered photons not in showers
    };

    struct Bounds
    {
        float xMin = std::numeric_limits<float>::max();
        float xMax = std::numeric_limits<float>::lowest();
        float zMin = std::numeric_limits<float>::max();
        float zMax = std::numeric_limits<float>::lowest();
        float rMax = 0;
        float xRange = 1.f;
        float zRange = 1.f;
    };

    pandora::StatusCode Run();

    /**
     *  @brief  Prepare training sample for the network
     *
     *  @return The StatusCode indicating success or failure in producing a training sample
     */
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief  Run network inference
     *
     *  @return The StatusCode indicating success or failure in running inference
     */
    pandora::StatusCode Infer();

    /**
     *  @brief  Fold EM showers to their leading particle.
     *
     *  @param mcHitMap The map of MC particles to their hits
     *  @param mcFoldingMap The map of MC particles to their hits
     *  @param leadingHitMap The map of leading MC particles to their consolidated hits
     */
    void FoldToLeading(const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, const MCFoldingMap &mcFoldingMap,
        lar_content::LArMCParticleHelper::MCContributionMap &leadingHitMap) const;

    /**
     *  @brief  Consolidate instances of hits from MC particles where folding is appropriate.
     *          Delta rays often own hits along muon or pion tracks, and these should be
     *          owned by the parent particle instead.
     *
     *  @param mcHitMap The map of MC particles to their hits
     *  @param instanceHitMap The map of MC particles to their consolidated hits
     */
    void ConsolidateInstances(const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap,
        lar_content::LArMCParticleHelper::MCContributionMap &instanceHitMap) const;

    /**
     *  @brief  Determine if an MC particle is a delta ray
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is a delta ray, false otherwise
     */
    bool IsDelta(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle is diffuse
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is diffuse, false otherwise
     */
    bool IsDiffuse(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle is highly ionising
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is highly ionising, false otherwise
     */
    bool IsHip(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle is minimally ionising
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is minimally ionising, false otherwise
     */
    bool IsMip(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle is a Michel electron
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is a Michel electron, false otherwise
     */
    bool IsMichel(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Determine if an MC particle is a shower
     *
     *  @param pLArMC The MC particle to check
     *  @return true if the MC particle is a shower, false otherwise
     */
    bool IsShower(const lar_content::LArMCParticle *const pLArMC) const;

    /**
     *  @brief  Choose the owner of a child MC particle based on its parent and the folding map
     *
     *  @param pMC The MC particle
     *  @param mcFoldingMap The map of MC particles to their folded counterparts
     *  @param mcHitMap The map of MC particles to their hits
     *  @param particleOwnedHits The list of hits owned by the particle
     *  @param parentOwnedHits The list of hits owned by the parent particle
     *
     *  @return The parent MC particle if there is allocation to a parent, nullptr otherwise
     */
    const lar_content::LArMCParticle *AllocateHitOwner(const pandora::MCParticle *const pMC, const MCFoldingMap &mcFoldingMap,
        const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, pandora::CaloHitList &particleOwnedHits,
        pandora::CaloHitList &parentOwnedHits) const;

    /**
     *  @brief  Choose the owner of a child MC particle based on its parent and the folding map
     *
     *  @param pMC The MC particle
     *  @param mcHitMap The map of MC particles to their hits
     *  @param particleOwnedHits The list of hits owned by the particle
     *  @param parentOwnedHits The list of hits owned by the parent particle
     *
     *  @return The parent MC particle if there is allocation to a parent, nullptr otherwise
     */
    const lar_content::LArMCParticle *AllocateHitOwner(const pandora::MCParticle *const pMC,
        const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, pandora::CaloHitList &particleOwnedHits,
        pandora::CaloHitList &parentOwnedHits) const;

    /**
     *  @brief  Get planar coordinates of vertex in given view
     *
     *  @param  view The view into which the vertex should be projected
     *  @param  vx The x position of the vertex
     *  @param  vz The z position of the vertex
     *
     *  @return The status code
     */
    pandora::StatusCode GetVertexPlanarCoordinates(pandora::HitType view, float &vx, float &vz) const;

    /**
     *  @brief  Get extrema of hit coordinates in the specified hit list
     *
     *  @param  caloHitList The list from which extrema should be extracted
     *  @param  vx The x position of the vertex
     *  @param  vz The z position of the vertex
     *  @param  bounds The output structure of extremal coordinate information
     */
    void GetCoordinateExtrema(const pandora::CaloHitList &caloHitList, [[maybe_unused]]const float vx, [[maybe_unused]]const float vz,
        Bounds &bounds) const;

    /**
     *  @brief  Populate input vectors for the deep learning network from the specified hit list
     *
     *  @param  caloHitList The list from which to populate the input vectors
     *  @param  bounds The coordinate bounds structure
     *  @param  vx The x position of the vertex
     *  @param  vz The z position of the vertex
     *  @param  xx The x coordinate vector
     *  @param  zz The z coordinate vector
     *  @param  rr The radial coordinate vector
     *  @param  cosTheta The cosine of the theta angle vector
     *  @param  sinTheta The sine of the theta angle vector
     *  @param  vv The value vector
     *  @param  adc The ADC value vector
     *  @param  width The width vector
     */
    void PopulateInputVectors(const pandora::CaloHitList &caloHitList, const Bounds &bounds, const float vx, const float vz, pandora::FloatVector &xx,
        pandora::FloatVector &zz, pandora::FloatVector &rr, pandora::FloatVector &cosTheta, pandora::FloatVector &sinTheta, pandora::FloatVector &vv,
        pandora::FloatVector &adc, pandora::FloatVector &width) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< Name of input calo hit list
    bool m_trainingMode; ///< Training mode
    bool m_vertexRelative; ///< Whether to make hit positions relative to the reconstructed vertex position
    bool m_polarCoords; ///< Whether to use polar coordinates for the hit positions
    float m_adcPeak; ///< Value representing peak of ADC distribution for all hits
    float m_maxAdcFactor; ///< Maximum ADC value for a hit will be m_adcPeak * m_maxAdcFactor, clipped beyond this
    int m_maxSeqLen; ///< Maximum sequence length when using polar coordinates
    std::string m_rootFileName; ///< Name of the ROOT file to save the training sample
    std::string m_rootTreeName; ///< Name of the ROOT tree to save the training sample
    std::string m_vertexListName; ///< Name of the vertex list to use for vertex-relative coordinates
    LArDLHelper::TorchModel m_model; ///< The model
};

} // namespace lar_dl_content

#endif // LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
