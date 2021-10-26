/**
 *  @file   larpandoracontent/LArMonitoring/AtmosphericsMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/AtmosphericsMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

AtmosphericsMonitoringAlgorithm::AtmosphericsMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AtmosphericsMonitoringAlgorithm::~AtmosphericsMonitoringAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_outputTreeName, m_outputFileName, "RECREATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "AtmosphericMonitoringAlgorithm: Unable to write to ROOT tree" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AtmosphericsMonitoringAlgorithm::Run()
{
    this->ConstructRootTree();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AtmosphericsMonitoringAlgorithm::ConstructRootTree() const
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
    const ParticleFlowObject *pRecoNeutrino{nullptr};
    if (!primaries.empty())
    {
        for (const MCParticle *primary : primaries)
        {
            const MCParticleList &parents{primary->GetParentList()};
            if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
            {
                pTrueNeutrino = parents.front();
                break;
            }
        }
    }

    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            pRecoNeutrino = pPfo;
            break;
        }
    }
    if (pTrueNeutrino && pRecoNeutrino)
    {
        const int flavour{pTrueNeutrino->GetParticleId()};
        const float energy{pTrueNeutrino->GetEnergy()};
        const int nhits{static_cast<int>(pCaloHitList->size())};

        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
        const CartesianVector &trueDirection{pTrueNeutrino->GetMomentum().GetUnitVector()};
        const CartesianVector tu(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ())));
        const CartesianVector tv(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ())));
        const CartesianVector tw(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ())));
        const float dirX{trueDirection.GetX()}, dirY{trueDirection.GetY()}, dirZ{trueDirection.GetZ()};

        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pRecoNeutrino)->GetPosition()};
        const CartesianVector ru(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(recoVertex.GetY(), recoVertex.GetZ())));
        const CartesianVector rv(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(recoVertex.GetY(), recoVertex.GetZ())));
        const CartesianVector rw(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(recoVertex.GetY(), recoVertex.GetZ())));

        const float dr{(ru - tu).GetMagnitude()};

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "flavour", flavour));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "energy", energy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "nhits", nhits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "vtx_dr", dr));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "dir_x", dirX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "dir_y", dirY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_outputTreeName, "dir_z", dirZ));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_outputTreeName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AtmosphericsMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFileName", m_outputFileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTreeName", m_outputTreeName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
