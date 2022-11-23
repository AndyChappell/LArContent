/**
 *  @file   larpandoracontent/LArMonitoring/CharacterisationMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of characterisation monitoring algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArMonitoring/CharacterisationMonitoringAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

CharacterisationMonitoringAlgorithm::CharacterisationMonitoringAlgorithm() :
    m_minHitsForGoodView{5}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CharacterisationMonitoringAlgorithm::~CharacterisationMonitoringAlgorithm()
{
    if (!m_rootTreeName.empty())
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "CharacterisationMonitoringAlgorithm: Unable to write to ROOT tree \'" << m_rootTreeName << "\' to file \'" << m_rootFileName <<
                "\'" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CharacterisationMonitoringAlgorithm::Run()
{
    for (const std::string &pfoListName : {m_trackPfoListName, m_showerPfoListName})
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            CaloHitList caloHitListU, caloHitListV, caloHitListW;
            const HitType viewU{HitType::TPC_VIEW_U}, viewV{HitType::TPC_VIEW_V}, viewW{HitType::TPC_VIEW_W};
            LArPfoHelper::GetCaloHits(pPfo, viewU, caloHitListU);
            LArPfoHelper::GetCaloHits(pPfo, viewV, caloHitListV);
            LArPfoHelper::GetCaloHits(pPfo, viewW, caloHitListW);
            int numGoodViews{0}, numHitsTotal{0};
            for (const CaloHitList &caloHitList : {caloHitListU, caloHitListV, caloHitListW})
            {
                numHitsTotal += caloHitList.size();
                if (caloHitList.size() >= 3)
                    ++numGoodViews;
            }
            if ((numGoodViews < 2) || (numHitsTotal < (3 * m_minHitsForGoodView)))
                continue;
            float showerContribution{0.f}, trackContribution{0.f};
            for (const CaloHitList &caloHitList : {caloHitListU, caloHitListV, caloHitListW})
            {
                for (const CaloHit *pCaloHit : caloHitList)
                {
                    try
                    {
                        const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                        if (!pMCParticle)
                            continue;
                        const int pdg{pMCParticle->GetParticleId()};
                        const int absPdg{std::abs(pdg)};
                        const float energy{pCaloHit->GetInputEnergy()};
                        if (absPdg == E_MINUS || absPdg == PHOTON)
                            showerContribution += energy;
                        else
                            trackContribution += energy;
                    }
                    catch (const StatusCodeException &)
                    {
                    }
                }
            }
            const float totalContribution{trackContribution + showerContribution};
            if (totalContribution <= std::numeric_limits<float>::epsilon())
                continue;
            const float trackFraction{trackContribution / totalContribution};
            const int isTrueTrack{trackFraction > 0.5f ? 1 : 0};
            const int isRecoTrack{pPfo->GetParticleId() == MU_MINUS};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "nHits", numHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "trueTrackFraction", trackFraction));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "isTrueTrack", isTrueTrack));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "isRecoTrack", isRecoTrack));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CharacterisationMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
    if (!m_rootFileName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

