/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/MyTestAlgorithm.cc
 *
 *  @brief  Implementation of the MyTest algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/MyTestAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MyTestAlgorithm::MyTestAlgorithm() :
    m_myMandatoryString(),
    m_myOptionalBool(false),
    m_myOptionalUnsignedInt(5),
    m_myMandatoryFloatVector()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyTestAlgorithm::Run()
{
    std::cout << "-m_myMandatoryString: " << m_myMandatoryString << std::endl
        << "-m_myOptionalBool: " << m_myOptionalBool << std::endl
        << "-m_myOptionalUnsignedInt: " << m_myOptionalUnsignedInt << std::endl
        << "-m_myMandatoryFloatVector: ";

    for (const auto value : m_myMandatoryFloatVector)
        std::cout << value << " ";

    std::cout << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyTestAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MyMandatoryString", m_myMandatoryString));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MyOptionalBool", m_myOptionalBool));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MyOptionalUnsignedInt", m_myOptionalUnsignedInt));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "MyMandatoryFloatVector", m_myMandatoryFloatVector));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

