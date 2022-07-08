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
     *  @brief  Add calo hits to the input for the network
     *
     *  @param  caloHitList The calo hits to be added
     *  @param  input The network input tensor to be populated
     */
    void PrepareNetworkInput(const pandora::CaloHitList &caloHitList, LArDLHelper::TorchInput &input) const;

    /**
     *  @brief  Identify the XZ range containing the hits for an event
     *
     *  @param  caloHitList The list of CaloHits for which the range is to be found
     *  @param  xMin The output minimum x-coordinate
     *  @param  xMax The output maximum x-coordinate
     *  @param  zMin The output minimum z-coordinate
     *  @param  zMax The output maximum z-coordinate
     */
    void GetHitRegion(const pandora::CaloHitList &caloHitList, double &xMin, double &xMax, double &zMin, double &zMax) const;

    /**
     *  @brief  Process a single PFO list
     */
    pandora::StatusCode ProcessPfoList(const std::string &pfoListName) const;

    /**
     *  @brief  Change the existing track/shower characterisation of a PFO
     *
     *  @param  pfo The particle flow object whose characterisation should be changed
     */
    pandora::StatusCode ChangeCharacterisation(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Adds hit information from a hit list to a feature vector
     */
    void PopulateFeatureVector(const pandora::CaloHitList &caloHitList, lar_content::LArMvaHelper::MvaFeatureVector &featureVector) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_training;    ///< Whether or not to run in training mode
    int m_imageWidth; ///< The width of event images
    int m_imageHeight; ///< The height of event images
    std::string m_trackPfoListName; ///< Name of input track PFO list
    std::string m_showerPfoListName; ///< Name of input shower PFO list
    std::string m_trainingFileName; ///< Name of the output training file
    std::string m_modelFileName; ///< Model file name
    LArDLHelper::TorchModel m_model; ///< Model
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H
