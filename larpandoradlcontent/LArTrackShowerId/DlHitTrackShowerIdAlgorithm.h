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

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< Name of input calo hit list
    bool m_trainingMode; ///< Training mode
    bool m_vertexRelative; ///< Whether to make hit positions relative to the reconstructed vertex position
    bool m_polarCoords; ///< Whether to use polar coordinates for the hit positions
    float m_adcNormalization; ///< ADC normalization factor - value representing peak of ADC distribution for all hits
    std::string m_rootFileName; ///< Name of the ROOT file to save the training sample
    std::string m_rootTreeName; ///< Name of the ROOT tree to save the training sample
    std::string m_vertexListName; ///< Name of the vertex list to use for vertex-relative coordinates
};

} // namespace lar_dl_content

#endif // LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
