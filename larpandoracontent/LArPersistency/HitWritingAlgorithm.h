/**
 *  @file   larpandoracontent/LArPersistency/HitWritingAlgorithm.h
 *
 *  @brief  Header file for the hit writing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_WRITING_ALGORITHM_H
#define LAR_HIT_WRITING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace pandora {class FileWriter;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  HitWritingAlgorithm class
 */
class HitWritingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HitWritingAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~HitWritingAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputCaloHitListNames;
    std::string m_inputCaloHitList2DName;
    std::string m_outputFilename;
    std::string m_treeName;

    const std::string m_interactions[31] = {
        "CCQEL_MU_P",
        "CCRES_MU_P_PIPLUS",
        "CCRES_MU_PIPLUS",
        "CCRES_MU_P_PIZERO",
        "CCRES_MU_P",
        "CCRES_MU",
        "CCRES_MU_P_P",
        "CCRES_MU_PIZERO",
        "CCRES_MU_P_P_P",
        "CCRES_MU_P_P_PIPLUS",
        "CCQEL_MU_P_P",
        "CCCOH",
        "CCRES_MU_P_P_P_P",
        "CCRES_MU_P_P_PIZERO",
        "CCQEL_MU_P_P_P",
        "CCQEL_E_P",
        "CCRES_MU_P_PHOTON",
        "CCRES_E_PIPLUS",
        "CCRES_MU_P_P_P_PIZERO",
        "CCRES_MU_P_P_P_P_P",
        "CCRES_E_P_PIPLUS",
        "CCRES_E_P",
        "CCQEL_MU_P_P_P_P",
        "CCRES_MU_PHOTON",
        "CCRES_MU_P_P_P_PIPLUS",
        "CCRES_MU_P_P_P_PHOTON",
        "CCRES_E_P_PHOTON",
        "CCRES_MU_P_P_P_P_P_PIZERO",
        "CCRES_E",
        "CCRES_MU_P_P_P_P_PIZERO",
        "CCRES_E_P_PIZERO"
    };
};

} // namespace lar_content

#endif // #ifndef LAR_HIT_WRITING_ALGORITHM_H
