/**
 *  @file   larpandoracontent/LArMonitoring/AdcMomentumAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArAtmos/AdcMomentumAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

AdcMomentumAlgorithm::AdcMomentumAlgorithm() :
    m_writeFile{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdcMomentumAlgorithm::~AdcMomentumAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AdcMomentumAlgorithm::Run()
{
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));


    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
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
    
    if (pTrueNeutrino && !pPfoList->empty())
    {
        std::cout << "Nu " << pTrueNeutrino->GetParticleId() << std::endl;
        const ParticleFlowObject *const pNeutrino{pPfoList->front()};
        this->MakePfoToMCMap(pNeutrino->GetDaughterPfoList());
           
        for (const ParticleFlowObject *const pPrimary : pNeutrino->GetDaughterPfoList())
        {
            const MCParticle *pMC{m_pfoToMCMap[pPrimary]};
            if (pMC->GetParentList().front()->GetParticleId() == 2112)
                continue;
            const double adcU{this->GetAdcContribution(pPrimary, TPC_VIEW_U)};
            const double adcV{this->GetAdcContribution(pPrimary, TPC_VIEW_V)};
            const double adcW{this->GetAdcContribution(pPrimary, TPC_VIEW_W)};
            const double adc{std::max({adcU, adcV, adcW})};
            const double momentum{pMC->GetMomentum().GetMagnitude()};
            const int mcPdg{pMC->GetParticleId()};

            if (m_writeFile)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mc_pdg", mcPdg));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mc_mom", momentum));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "pfo_adc", adc));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
            }
            std::cout << mcPdg << " " << momentum << " " << adc << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdcMomentumAlgorithm::GetAdcContribution(const ParticleFlowObject *const pPfo, const HitType view) const
{
    double totalAdc{0.};

    CaloHitList caloHits;
    LArPfoHelper::GetCaloHits(pPfo, view, caloHits);

    for(const CaloHit *const pCaloHit : caloHits)
        totalAdc += pCaloHit->GetInputEnergy();

    for (const ParticleFlowObject *const pChildPfo : pPfo->GetDaughterPfoList())
        totalAdc += this->GetAdcContribution(pChildPfo, view);

    return totalAdc;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AdcMomentumAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AdcMomentumAlgorithm::MakePfoToMCMap(const PfoList &pfoList)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        CaloHitList caloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, caloHits);
        std::map<const MCParticle*, int> m_mcToHitMap;
        for (const CaloHit *const pCaloHit : caloHits)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (m_mcToHitMap.find(pMC) != m_mcToHitMap.end())
                    ++m_mcToHitMap[pMC];
                else
                    m_mcToHitMap[pMC] = 1;
            }
            catch (...)
            {
            }
        }
        int bestCount{0};
        m_pfoToMCMap[pPfo] = nullptr;
        for (const auto & [pMC, count] : m_mcToHitMap)
        {
            if (count > bestCount)
            {
                bestCount = count;
                m_pfoToMCMap[pPfo] = pMC;
            }
        }
    }
}

} // namespace lar_content
