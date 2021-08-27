/**
 *  @file   larpandoracontent/LArControlFlow/HitWritingAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/HitWritingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

HitWritingAlgorithm::HitWritingAlgorithm()
{
}

HitWritingAlgorithm::~HitWritingAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "hits", "hits.root", "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWritingAlgorithm::Run()
{
    const CaloHitList *pCaloHitListU(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListU", pCaloHitListU));
    const CaloHitList *pCaloHitListV(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListV", pCaloHitListV));
    const CaloHitList *pCaloHitListW(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitListW));

    FloatVector xHitsU, zHitsU, xHitsV, zHitsV, xHitsW, zHitsW;
    for (const CaloHit *pCaloHit : *pCaloHitListU)
    {
        xHitsU.emplace_back(pCaloHit->GetPositionVector().GetX());
        zHitsU.emplace_back(pCaloHit->GetPositionVector().GetZ());
    }
    for (const CaloHit *pCaloHit : *pCaloHitListV)
    {
        xHitsV.emplace_back(pCaloHit->GetPositionVector().GetX());
        zHitsV.emplace_back(pCaloHit->GetPositionVector().GetZ());
    }
    for (const CaloHit *pCaloHit : *pCaloHitListW)
    {
        xHitsW.emplace_back(pCaloHit->GetPositionVector().GetX());
        zHitsW.emplace_back(pCaloHit->GetPositionVector().GetZ());
    }
    std::cout << "List sizes: " << xHitsU.size() << " " << xHitsV.size() << " " << xHitsW.size() << std::endl;
  
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "hits", "x_u", &xHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "hits", "z_u", &zHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "hits", "x_v", &xHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "hits", "z_v", &zHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "hits", "x_w", &xHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "hits", "z_w", &zHitsW));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "hits"));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWritingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
