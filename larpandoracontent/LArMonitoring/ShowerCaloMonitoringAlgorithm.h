/**
 *  @file   larpandoracontent/LArMonitoring/ShowerCaloMonitoringAlgorithm.h
 *
 *  @brief  Header file for the shower calo algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_CALO_MONITORING_ALGORITHM_H
#define LAR_SHOWER_CALO_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <memory>

namespace lar_content
{

/**
 *  @brief  ShowerCaloMonitoringAlgorithm class
 */
class ShowerCaloMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerCaloMonitoringAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~ShowerCaloMonitoringAlgorithm();

private:
    typedef std::pair<const float, const pandora::CaloHit *> ProjCaloHit;
    typedef std::shared_ptr<ProjCaloHit> ProjCaloHitPtr;
    typedef std::vector<ProjCaloHitPtr> ProjCaloHitPtrList;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MakeSelection(const pandora::MCParticleList *pMCList, const pandora::CaloHitList *pCaloHitList,
        LArMCParticleHelper::MCContributionMap &mcMap) const;

    void ProjectHit(const pandora::CaloHit *pCaloHit, const pandora::CartesianVector &origin, const pandora::CartesianVector &axis,
        ProjCaloHitPtrList &projCaloHitList);

    std::string m_caloHitListName; ///< Name of input calo hit list
    bool m_visualize;              ///< Whether or not to visualize MC particles
};

} // namespace lar_content

#endif // LAR_SHOWER_CALO_MONITORING_ALGORITHM_H
