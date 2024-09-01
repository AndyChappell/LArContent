/**
 *  @file   larpandoracontent/LArMonitoring/MCHierarchyAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_HIERARCHY_ALGORITHM_H
#define LAR_MC_HIERARCHY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  MCHierarchyAlgorithm class
 */
class MCHierarchyAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCHierarchyAlgorithm();

    ~MCHierarchyAlgorithm();

private:
    pandora::StatusCode Run();

    void Make3DHits(const pandora::CaloHitList &hits2D, pandora::CaloHitList &hits3D) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< Name of input calo hit list
    std::string m_mcParticleListName;  ///< Name of input MC particle list
    std::string m_mcFileName; ///< Name of the MC output file
    std::string m_mcTreeName; ///< Name of the MC output tree
    std::string m_eventFileName; ///< Name of the event output file
    std::string m_eventTreeName; ///< Name of the event output tree
    std::string m_hitsFileName; ///< Name of the hits output file
    std::string m_hitsTreeName; ///< Name of the hits output tree
    float m_correctionX; ///< Truth correction to apply in x
    float m_correctionY; ///< Truth correction to apply in y
    float m_correctionZ; ///< Truth correction to apply in z
    bool m_visualize; // Whether to visualize the event
    LArMCParticleHelper::MCContributionMap m_mcToHitsMap;
};

} // namespace lar_content

#endif // LAR_MC_HIERARCHY_ALGORITHM_H
