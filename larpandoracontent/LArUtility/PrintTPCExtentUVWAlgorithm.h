/**
 *  @file   larpandoracontent/LArUtility/PrintTPCExtentAlgorithm.h
 *
 *  @brief  Header file for the TPC extent printing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TPC_EXTENT_UVW_ALGORITHM_H
#define LAR_TPC_EXTENT_UVW_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PrintTPCExtentAlgorithm class
 */
class PrintTPCExtentUVWAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_TPC_EXTENT_UVW_ALGORITHM_H
