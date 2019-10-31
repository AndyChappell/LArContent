/**
 *  @file   larpandoracontent/LArUtility/PrintTPCExtentUVWAlgorithm.cc
 *
 *  @brief  Implementation of the list merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/PrintTPCExtentUVWAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode PrintTPCExtentUVWAlgorithm::Run()
{
    const Pandora& pandora = this->GetPandora();
    const LArTPCMap& larTPCMap(pandora.GetGeometry()->GetLArTPCMap());
    
    // Currently assumes a single TPC in the map
    float xMin(0.f), xMax(0.f), yMin(0.f), yMax(0.f), zMin(0.f), zMax(0.f);
    for(auto iter : larTPCMap)
    {
        const LArTPC* const larTPC(iter.second);
        const float xWidth2 = larTPC->GetWidthX() * 0.5f;
        xMin = larTPC->GetCenterX() - xWidth2;
        xMax = larTPC->GetCenterX() + xWidth2;
        const float yWidth2 = larTPC->GetWidthY() * 0.5f;
        yMin = larTPC->GetCenterY() - yWidth2;
        yMax = larTPC->GetCenterY() + yWidth2;
        const float zWidth2 = larTPC->GetWidthZ() * 0.5f;
        zMin = larTPC->GetCenterZ() - zWidth2;
        zMax = larTPC->GetCenterZ() + zWidth2;
    }
    
    const float uMin = pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(yMin, zMin);
    const float uMax = pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoU(yMax, zMax);
    const float vMin = pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(yMin, zMin);
    const float vMax = pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoV(yMax, zMax);
    const float wMin = pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(yMin, zMin);
    const float wMax = pandora.GetPlugins()->GetLArTransformationPlugin()->YZtoW(yMax, zMax);
    
    std::cout << "TPC min extent XYZ(" << xMin << "," << yMin << "," << zMin << ") => XUVW(" <<
        xMin << "," << uMin << "," << vMin << "," << wMin << ")" << std::endl;
    std::cout << "TPC max extent XYZ(" << xMax << "," << yMax << "," << zMax << ") => XUVW(" <<
        xMax << "," << uMax << "," << vMax << "," << wMax << ")" << std::endl;    
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PrintTPCExtentUVWAlgorithm::ReadSettings(const TiXmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
