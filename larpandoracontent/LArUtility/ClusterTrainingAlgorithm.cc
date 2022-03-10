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
            CaloHitList volumeHitsU, volumeHitsV, volumeHitsW;
            volume.GetHitList(TPC_VIEW_U, volumeHitsU);
            volume.GetHitList(TPC_VIEW_V, volumeHitsV);
            volume.GetHitList(TPC_VIEW_W, volumeHitsW);
            this->CreateHitList(volumeHitsU, "CaloHitListU_0_" + std::to_string(t) + "_" + std::to_string(c));
            this->CreateHitList(volumeHitsV, "CaloHitListV_0_" + std::to_string(t) + "_" + std::to_string(c));
            this->CreateHitList(volumeHitsW, "CaloHitListW_0_" + std::to_string(t) + "_" + std::to_string(c));
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

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterTrainingAlgorithm::CreateHitList(const CaloHitList &caloHitList, std::string listName) const
{
    return PandoraContentApi::SaveList<CaloHitList>(*this, caloHitList, listName);
}

} // namespace lar_content

