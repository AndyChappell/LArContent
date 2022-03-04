/**
 *  @file   larpandoracontent/LArUtility/ClusterTrainingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster training algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArTpcGeometryHelper.h"
#include "larpandoracontent/LArUtility/ClusterTrainingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterTrainingAlgorithm::ClusterTrainingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterTrainingAlgorithm::Run()
{
    const LArTransformationPlugin *const pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    LArTpcGeometryHelper helper{LArTpcGeometryHelper::GetInstance(pTransform)};
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));
    for (unsigned int t = 0; t < 2; ++t)
    {
        for (unsigned int c = t % 2; c < 24; c += 2)
        {
            LArTpcGeometryHelper::VolumeId id(0, t, c);
            TpcHitVolume &volume{helper.GetTpcHitVolume(id)};
            volume.Add(*pCaloHitList);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterTrainingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    (void)xmlHandle;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

