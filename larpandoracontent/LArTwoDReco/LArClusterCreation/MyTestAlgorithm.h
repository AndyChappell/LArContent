/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/MyTestAlgorithm.h
 *
 *  @brief  Header file for the MyTest algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MY_TEST_ALGORITHM_H
#define LAR_MY_TEST_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MyTestAlgorithm class
 */
class MyTestAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MyTestAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
    std::string m_myMandatoryString; ///< A mandatory string
    bool m_myOptionalBool; ///< An optional boolean
    unsigned int m_myOptionalUnsignedInt; ///< An optional unsigned int
    pandora::FloatVector m_myMandatoryFloatVector; ///< A mandatory vector of floats
};

} // namespace lar_content

#endif // #ifndef LAR_MY_TEST_ALGORITHM_H
