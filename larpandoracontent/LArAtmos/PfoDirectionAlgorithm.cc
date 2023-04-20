/**
 *  @file   larpandoracontent/LArMonitoring/PfoDirectionAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArAtmos/PfoDirectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

PfoDirectionAlgorithm::PfoDirectionAlgorithm() :
    m_writeFile{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoDirectionAlgorithm::~PfoDirectionAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoDirectionAlgorithm::Run()
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
    MCParticleList primariesList(primaries.begin(), primaries.end());
    const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primariesList)};
    const int isCC{descriptor.IsCC()};
    const int isQE{descriptor.IsQE()};
    const int isRes{descriptor.IsResonant()};
    const int isDIS{descriptor.IsDIS()};
    const int isCoh{descriptor.IsCoherent()};
    const int isMu{descriptor.IsMuonNeutrino()};
    const int isElectron{descriptor.IsElectronNeutrino()};

    CartesianVector mcDirection(0, 0, 0);
    for (const MCParticle *const pMC : primaries)
    {
        std::cout << "MC " << pMC->GetParticleId() << " " << pMC->GetMomentum().GetUnitVector() << std::endl;
        mcDirection += pMC->GetMomentum();
    }
    mcDirection = mcDirection.GetUnitVector();
    std::cout << "MC Direction: " << mcDirection << std::endl;

    if (pTrueNeutrino && !pPfoList->empty())
    {
        std::map<const ParticleFlowObject *const, double> pfoAdcContribution;
        double totalAdc{0.};
        const ParticleFlowObject *const pNeutrino{pPfoList->front()};
        this->MakePfoToMCMap(pNeutrino->GetDaughterPfoList());
           
        for (const ParticleFlowObject *const pPrimary : pNeutrino->GetDaughterPfoList())
        {
            CaloHitList caloHits;
            LArPfoHelper::GetCaloHits(pPrimary, TPC_3D, caloHits);
            const MCParticle *pMC{m_pfoToMCMap[pPrimary]};
            if (caloHits.empty() || (pMC->GetParentList().front()->GetParticleId() == 2112))
                continue;
            const double adcU{this->GetAdcContribution(pPrimary, TPC_VIEW_U)};
            const double adcV{this->GetAdcContribution(pPrimary, TPC_VIEW_V)};
            const double adcW{this->GetAdcContribution(pPrimary, TPC_VIEW_W)};
            pfoAdcContribution[pPrimary] = std::max({adcU, adcV, adcW});
            //pfoAdcContribution[pPrimary] = std::sqrt(this->GetAdcContribution(pPrimary));
            totalAdc += pfoAdcContribution[pPrimary];

            // Get distribution for momentum magnitude against adc for different particle types to see if there is
            // some useful correction that can be applied
        
            std::cout << "PFO " << pPrimary->GetParticleId() << " Hits: " << caloHits.size() << " (" << pfoAdcContribution[pPrimary] << ") " <<
                pPrimary->GetMomentum().GetUnitVector() << std::endl;
            std::cout << "   MC " << pMC->GetParticleId() << " " << pMC->GetMomentum() << std::endl;
            std::cout << "   " << (pfoAdcContribution[pPrimary] / pMC->GetMomentum().GetMagnitude()) << " " <<
                (pfoAdcContribution[pPrimary] / pMC->GetMomentum().GetMagnitudeSquared()) << std::endl;
        }

        CartesianVector recoDirection(0, 0, 0);
        for (const auto & [ pPfo, adc ] : pfoAdcContribution)
            recoDirection += pPfo->GetMomentum().GetUnitVector() * (adc / totalAdc);
        recoDirection = recoDirection.GetUnitVector();

        const CartesianVector &trueDirection{pTrueNeutrino->GetMomentum().GetUnitVector()};
        std::cout << "True energy: " << pTrueNeutrino->GetEnergy() << std::endl;
        std::cout << (isCC ? "CC" : "NC") << (isQE ? "QE" : (isRes ? "RES" : (isDIS ? "DIS" : "Other"))) << std::endl;
        std::cout << "True Direction: " << trueDirection << std::endl;
        std::cout << "Reco Direction: " << recoDirection << std::endl;

        const int truePdg{pTrueNeutrino->GetParticleId()};
        if (m_writeFile)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_pdg", truePdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_dir_x", trueDirection.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_dir_y", trueDirection.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "true_dir_z", trueDirection.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mc_dir_x", mcDirection.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mc_dir_y", mcDirection.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mc_dir_z", mcDirection.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "reco_dir_x", recoDirection.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "reco_dir_y", recoDirection.GetY()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "reco_dir_z", recoDirection.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "is_cc", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "is_qe", isQE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "is_res", isRes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "is_dis", isDIS));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "is_coh", isCoh));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "is_mu", isMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "is_e", isElectron));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double PfoDirectionAlgorithm::GetAdcContribution(const ParticleFlowObject *const pPfo, const HitType view) const
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

StatusCode PfoDirectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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

void PfoDirectionAlgorithm::MakePfoToMCMap(const PfoList &pfoList)
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
