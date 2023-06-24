/**
 *  @file   larpandoracontent/LArWorkshop/MyClusteringAlgorithm.cc
 *
 *  @brief  Implementation of a custom clustering algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArWorkshop/MyClusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MyClusteringAlgorithm::MyClusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusteringAlgorithm::Run()
{
    // Your algorithm goes here

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    (void)xmlHandle;
    // Read in any XML parameters here

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

