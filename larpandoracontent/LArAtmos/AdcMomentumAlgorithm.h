/**
 *  @file   larpandoracontent/LArMonitoring/AdcMomentumAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_ADC_MOMENTUM_ALGORITHM_H
#define LAR_ADC_MOMENTUM_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  AdcMomentumAlgorithm class
 */
class AdcMomentumAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    AdcMomentumAlgorithm();

    virtual ~AdcMomentumAlgorithm();

private:
    pandora::StatusCode Run();
    double GetAdcContribution(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType view) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MakePfoToMCMap(const pandora::PfoList &pfoList);

    bool m_writeFile;               // Whether to produce ROOT output file
    std::string m_filename;         // The filename of the ROOT output file
    std::string m_treename;         // The name of the ROOT tree

    std::map<const pandora::ParticleFlowObject*, const pandora::MCParticle*> m_pfoToMCMap;
};

} // namespace lar_content

#endif // LAR_ADC_MOMENTUM_ALGORITHM_H
