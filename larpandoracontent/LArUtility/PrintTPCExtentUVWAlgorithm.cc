/**
 *  @file   larpandoracontent/LArUtility/PrintTPCExtentUVWAlgorithm.cc
 *
 *  @brief  Implementation of the list merging algorithm class.
 *
 *  $Log: $
 */

#include <algorithm>
#include <limits>
#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/PandoraEnumeratedTypes.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

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
        xMin = std::min(xMin, larTPC->GetCenterX() - xWidth2);
        xMax = std::max(xMax, larTPC->GetCenterX() + xWidth2);
        const float yWidth2 = larTPC->GetWidthY() * 0.5f;
        yMin = std::min(yMin, larTPC->GetCenterY() - yWidth2);
        yMax = std::max(yMax, larTPC->GetCenterY() + yWidth2);
        const float zWidth2 = larTPC->GetWidthZ() * 0.5f;
        zMin = std::min(zMin, larTPC->GetCenterZ() - zWidth2);
        zMax = std::max(zMax, larTPC->GetCenterZ() + zWidth2);
        
        std::cout << "X: " << larTPC->GetCenterX() << " " << larTPC->GetWidthX() << std::endl;
        std::cout << "Y: " << larTPC->GetCenterY() << " " << larTPC->GetWidthY() << std::endl;
        std::cout << "Z: " << larTPC->GetCenterZ() << " " << larTPC->GetWidthZ() << std::endl;
    }
    
    const CartesianVector min(xMin, yMin, zMin);
    const CartesianVector max(xMax, yMax, zMax);
    
    float uMin(std::numeric_limits<float>::max()), uMax(-std::numeric_limits<float>::max());
    float vMin(uMin), vMax(uMax), wMin(uMin), wMax(uMax);
    auto plugin = pandora.GetPlugins()->GetLArTransformationPlugin();
    for(float y = min.GetY(); y < max.GetY() + 0.45f; y += 0.45f)
    {
        for(float z = min.GetZ(); z < max.GetZ() + 0.45f; z += 0.45f)
        {
            uMin = std::min(uMin, static_cast<float>(plugin->YZtoU(y, z)));
            uMax = std::max(uMax, static_cast<float>(plugin->YZtoU(y, z)));
            vMin = std::min(vMin, static_cast<float>(plugin->YZtoV(y, z)));
            vMax = std::max(vMax, static_cast<float>(plugin->YZtoV(y, z)));
            wMin = std::min(wMin, static_cast<float>(plugin->YZtoW(y, z)));
            wMax = std::max(wMax, static_cast<float>(plugin->YZtoW(y, z)));
        }
    }

    std::cout << "Min XUVW(" << xMin << "," << uMin << "," << vMin << "," << wMin << ")" << std::endl;
    std::cout << "Max XUVW(" << xMax << "," << uMax << "," << vMax << "," << wMax << ")" << std::endl;

    //CartesianVector uPosMin = LArGeometryHelper::ProjectPosition(pandora, posMin, TPC_VIEW_U);

    /*
    CartesianVector minTest(158.f, -102.653f, 173.038f);
    CartesianVector maxTest(1353.5f, 1451.28f, 1391.99f);
    
    const float YfromUV(pandora.GetPlugins()->GetLArTransformationPlugin()->UVtoY(uPosMin.GetZ(), vPosMin.GetZ()));
    const float YfromUW(pandora.GetPlugins()->GetLArTransformationPlugin()->UWtoY(uPosMin.GetZ(), wPosMin.GetZ()));
    const float YfromVW(pandora.GetPlugins()->GetLArTransformationPlugin()->VWtoY(vPosMin.GetZ(), wPosMin.GetZ()));
    
    std::cout << "Y from UV: " << YfromUV << "   from UW: " << YfromUW << "   from VW: " << YfromVW << std::endl;
    */
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PrintTPCExtentUVWAlgorithm::ReadSettings(const TiXmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
