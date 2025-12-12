/**
 *  @file   larpandoracontent/LArHelpers/LArTopologyHelper.h
 *
 *  @brief  Header file for the principal curve analysis helper class.
 *
 *  $Log: $
 */
#ifndef LAR_TOPOLOGY_HELPER_H
#define LAR_TOPOLOGY_HELPER_H 1

//#include "Api/PandoraApi.h"
#include "Objects/CaloHit.h"
#include "Objects/CartesianVector.h"
#include "Objects/MCParticle.h"
#include "Plugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <map>

namespace lar_content
{

/**
 *  @brief  LArTopologyHelper class
 */
class LArTopologyHelper
{
public:
    typedef std::map<const pandora::MCParticle*, pandora::CartesianPointVector> MCVertexMap;

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

    struct Filter
    {
        /**
         *  @brief  Constructor
         *
         *  @param  minHits the minimum number of hits required to retain a vertex (not applied to electrons and photons)
         *  @param  neutrons whether to filter out neutron-induced activity
         *  @param  deltas whether to filter out delta rays
         *  @param  compton whether to filter out Compton scattering
         *  @param  photoElectric whether to filter out photoelectric effect
         */
        Filter(const size_t minHits = 0, const bool neutrons = true, const bool deltas = true, const bool compton = true, const bool photoElectric = false);

        const size_t m_minHits;
        const bool m_neutrons;
        const bool m_deltas;
        const bool m_compton;
        const bool m_photoElectric;
    };

    typedef std::map<const pandora::CartesianVector, pandora::CaloHitList, VertexComparator> VertexHitsMap;
    typedef std::map<const pandora::CaloHit*, int> CondensationPointMap;
    typedef std::map<const pandora::CaloHit*, pandora::IntVector> HitVertexLabelMap;
    typedef std::map<const pandora::CaloHit*, pandora::FloatVector> HitVertexWeightMap;

    /**
     *  @brief  Get the topological vertices from a collection of hits. Here we define a topological vertex as a true vertex projected onto a TPC,
     *          view, potentially filtering out neutron-induced activity. The intent is for this to broadly correspond to the vertices you would
     *          identify by eye.
     *
     *  @param  caloHitList the collection of hits
     *  @param  filter the filtering options
     *  @param  pTransform the LAr transformation plugin
     *  @param  vertices the output collection of topological vertices
     */
    static void GetTopologicalVertices(const pandora::CaloHitList &caloHitList, const Filter &filter,
        const pandora::LArTransformationPlugin *const pTransform, pandora::CartesianPointVector &vertices);

    /**
     *  @brief  Get the projected true vertex positions in a given view
     *
     *  @param  pMC the MC particle for which the vertex should be found
     *  @param  view the TPC view
     *  @param  pTransform the LAr transformation plugin
     *  @param  mcVertexMap the output mapping between MC particles and projected true vertex positions
     */
    static void GetProjectedTrueVertices(const pandora::MCParticle *const pMC, const pandora::HitType view,
        const pandora::LArTransformationPlugin *const pTransform,MCVertexMap &mcVertexMap);

    /**
     *  @brief  Match each vertex to the closest hit and return the distance
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     *  @param  mcVertexMap the mapping between MC particles and their projected true vertex positions
     *  @param  mcToMatchedVertexMap the output mapping between MC particles and their matched vertex positions
     */
    static void MatchHitToVertex(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCVertexMap &mcVertexMap, MCVertexMap &mcToMatchedVertexMap);

    /**
     *  @brief  Consolidate any duplicate vertices and associate the relevant hits to the consolidated vertex
     *
     *  @param  mcToMatchedVertexMap the mapping between MC particles and their matched vertex positions
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     *  @param  vertexHitsMap the output mapping between consolidated vertex positions and their associated hits
     */
    static void ConsolidateVertices(const MCVertexMap &mcToMatchedVertexMap, const LArMCParticleHelper::MCContributionMap &mcToHitsMap, VertexHitsMap &vertexHitsMap);

    /**
     *  @brief  General filtering function
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     *  @param  filter the filtering options
     */
    static void FilterVertices(LArMCParticleHelper::MCContributionMap &mcToHitsMap, const Filter &filter);

    /**
     *  @brief  Filter out delta ray vertices
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     */
    static void FilterDeltaRays(LArMCParticleHelper::MCContributionMap &mcToHitsMap);

    /**
     *  @brief  Filter out particles that are neutron-induced
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     */
    static void FilterNeutronInducedParticles(LArMCParticleHelper::MCContributionMap &mcToHitsMap);

    /**
     *  @brief  Filter out Compton scattering vertices
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     */
    static void FilterComptonScattering(LArMCParticleHelper::MCContributionMap &mcToHitsMap);

    /**
     *  @brief  Filter out photoelectron vertices
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     */
    static void FilterPhotoElectrons(LArMCParticleHelper::MCContributionMap &mcToHitsMap);

    /**
     *  @brief  Filter out vertices with insufficient hits
     *
     *  @param  mcToHitsMap the mapping between MC particles and their contributed hits
     *  @param  minHits the minimum number of hits required to retain a vertex
     */
    static void FilterByMinHits(LArMCParticleHelper::MCContributionMap &mcToHitsMap, const size_t minHits);
};

} // namespace lar_content

#endif // #ifndef LAR_TOPOLOGY_HELPER_H
