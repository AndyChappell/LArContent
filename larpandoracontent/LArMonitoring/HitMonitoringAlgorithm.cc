/**
 *  @file   larpandoracontent/LArMonitoring/HitMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/HitMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

HitMonitoringAlgorithm::HitMonitoringAlgorithm() :
    m_caloHitListName{"CaloHitList2D"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitMonitoringAlgorithm::~HitMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename, m_filename, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitMonitoringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (pCaloHitList && !pCaloHitList->empty())
    {
        this->GetHitWidthDistribution(*pCaloHitList);
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitMonitoringAlgorithm::GetHitWidthDistribution(const CaloHitList &caloHitList) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const int view{pCaloHit->GetHitType()};
        const float width{pCaloHit->GetCellSize1()};
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "view", view));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename, "width", width));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_filename));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_treename));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

