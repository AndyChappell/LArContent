/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{

/**
 *  @brief  DlPfoCharacterisationAlgorithm class
 */
class DlPfoCharacterisationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlPfoCharacterisationAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Run the network inference to characterise the PFOs
     */
    pandora::StatusCode Infer();

    /**
     *  @brief  Produce files that act as inputs to network training
     */
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief  Prepares input for the network based on the provided PfoList
     *
     *  @param  pfoList The calo hits to be added
     *  @param  input The network input tensor to be populated
     */
    void PrepareNetworkInput(const pandora::PfoList &pfoList, LArDLHelper::TorchInput &input) const;

    /**
     *  @brief  Process a single PFO list
     */
    pandora::StatusCode ProcessPfoList(const std::string &pfoListName) const;

    /**
     *  @brief  Process a single PFO view
     *
     *  @param  caloHitList The CaloHitList from a single view
     *  @param  featureVector The feature vector to populate
     */
    void ProcessPfoView(const pandora::CaloHitList &caloHitList, lar_content::LArMvaHelper::MvaFeatureVector &featureVector) const;

    /**
     *  @brief Retrieves the fraction of hit energy in a PFO contributed from track-like hits
     *
     *  @param pPfo The PFO to consider
     *
     *  @returns The fraction of hit energy coming from track-like hits
     */
    float GetTrackFraction(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Change the existing track/shower characterisation of a PFO
     *
     *  @param  pfo The particle flow object whose characterisation should be changed
     */
    pandora::StatusCode ChangeCharacterisation(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Adds hit information from a hit list to a feature vector
     */
    void PopulateFeatureVector(const pandora::FloatVector &longitudinalProfile, const pandora::FloatVector &transverseProfile,
        lar_content::LArMvaHelper::MvaFeatureVector &featureVector) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_visualize;   ///< Wehter or not to visualize the charactierisation
    bool m_training;    ///< Whether or not to run in training mode
    std::string m_trackPfoListName; ///< Name of input track PFO list
    std::string m_showerPfoListName; ///< Name of input shower PFO list
    std::string m_trainingFileName; ///< Name of the output training file
    std::string m_modelFileName; ///< Model file name
    LArDLHelper::TorchModel m_model; ///< Model
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H
