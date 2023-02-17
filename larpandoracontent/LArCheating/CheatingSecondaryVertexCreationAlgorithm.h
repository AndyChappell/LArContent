/**
 *  @file   larpandoracontent/LArCheating/CheatingSecondaryVertexCreationAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SECONDARY_VERTEX_CREATION_ALGORITHM_H
#define LAR_CHEATING_SECONDARY_VERTEX_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingSecondaryVertexCreationAlgorithm::Algorithm class
 */
class CheatingSecondaryVertexCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingSecondaryVertexCreationAlgorithm();

private:
    typedef std::map<const pandora::MCParticle *, bool> VisibleParticleMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Identify the secondary vertex if it exists. For a secondary vertex to exist in this case we want a visible primary to interact
     *          and produce at least one visible secondary particle, where there is at least one secondary with a distinct PDG, or, if no
     *          distinct PDG code, the change in trajectory must be obvious by eye.
     *
     *  @param  pParent The parent MC particle to use as the starting point for identifying secondaries
     *  @param  visibleParticleMap A map indicating all visible MC particles
     *  @param  vertex The output vertex position
     *
     *  @return true if a visible secondary particle was found, false otherwise
     */
    bool GetClearSecondaryVertex(const pandora::MCParticle *const pParent, const VisibleParticleMap &visibleParticleMap,
        pandora::CartesianVector &vertex) const;

    std::string m_outputVertexListName; ///< The name under which to save the output vertex list
    bool m_replaceCurrentVertexList;    ///< Whether to replace the current vertex list with the output list
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SECONDARY_SECONDARY_VERTEX_CREATION_ALGORITHM_H
