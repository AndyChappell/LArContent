/**
 *  @file   larpandoradlcontent/LArThreeDReco/DlPfoEvaluatorAlgorithm.h
 *
 *  @brief  Header file for the deep learning PFO evaluator algorithm. Assesses the reconstruction quality of PFOs to identify candidates for
 *          reclustering.
 *
 *  $Log: $
 */
#ifndef LAR_DL_PFO_EVALUATOR_ALGORITHM_H
#define LAR_DL_PFO_EVALUATOR_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <map>

namespace lar_dl_content
{

/**
 *  @brief  DlPfoEvaluatorAlgorithm class
 */
class DlPfoEvaluatorAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlPfoEvaluatorAlgorithm();

    ~DlPfoEvaluatorAlgorithm();

private:
    typedef std::map<const pandora::MCParticle *, pandora::CaloHitList> MCHitMap;
    typedef std::map<pandora::HitType, MCHitMap> MapOfMCHitMaps;

    pandora::StatusCode Run();
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();
    float GetPurity(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType view);
    float GetCompleteness(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType view, const MapOfMCHitMaps &mcToAllHitsMap);
    float GetFactor(const pandora::CaloHitList &caloHitList);
    void CreateTrainingExample(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType view, const float metric, const float threshold);
    void GetMCHitsMap(const pandora::CaloHitList caloHitList, MCHitMap &mcHitsMap);
    void GetMCHitsMap(const pandora::CaloHitList caloHitList, MapOfMCHitMaps &mcAllHitsMap);
    const pandora::MCParticle *GetMainMCParticle(const pandora::CaloHitList &pfoHitList, float &weight);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;  ///< The names of the PFO lists to evaluate
    std::string m_caloHitListName;  ///< The name of the CaloHitList containing all 2D hits
    std::string m_rootFileName; ///< The output ROOT filename
    std::string m_rootTreePrefix; ///< The prefix for the ROOT tree name
    float m_purityThreshold;    ///< The threshold above which a PFO is considered well-reconstructed based on purity
    float m_completenessThreshold;    ///< The threshold above which a PFO is considered well-reconstructed based on completeness
    bool m_useAdcWeighting; ///< Whether or not hit contributions should be weighted by ADC
    bool m_trainingMode;    ///< Whether or not the algorithm is in training mode
};

} // namespace lar_dl_content

#endif // LAR_DL_PFO_EVALUATOR_ALGORITHM_H
