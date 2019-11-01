/**
 *  @file   larpandoracontent/LArPersistency/BScProjectInputsAlgorithm.h
 *
 *  @brief  Header file for the BSc project inputs algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_BSC_PROJECT_INPUTS_ALGORITHM_H
#define LAR_BSC_PROJECT_INPUTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <fstream>

namespace pandora {class FileWriter;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  BScProjectInputsAlgorithm class
 */
class BScProjectInputsAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    BScProjectInputsAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~BScProjectInputsAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    std::string m_inputCaloHitListName;
    std::string m_inputCaloHitListNameW;
    std::string m_outputFilename;
    std::ofstream m_file;
    
};

} // namespace lar_content

#endif // #ifndef LAR_BSC_PROJECT_INPUTS_ALGORITHM_H
