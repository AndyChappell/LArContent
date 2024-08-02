/**
 *  @file   larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

VertexMonitoringAlgorithm::VertexMonitoringAlgorithm() :
    m_event{-1},
    m_visualise{true},
    m_writeFile{false},
    m_fromList{false},
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexMonitoringAlgorithm::~VertexMonitoringAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::Run()
{
    ++m_event;
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(
            this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));
    }

    this->AssessVertices();

    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexMonitoringAlgorithm::Print(const MCParticle *const pMC, const std::string &indent) const
{
    if (pMC)
    {
        std::cout << indent << pMC->GetParticleId() << ": " << pMC->GetEnergy() << std::endl;
        for (const MCParticle *const pChild : pMC->GetDaughterList())
        {
            this->Print(pChild, indent + "  ");
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::AssessVertices() const
{
#ifdef MONITORING
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    int nHitsU{0}, nHitsV{0}, nHitsW{0};
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                ++nHitsU;
                break;
            case TPC_VIEW_V:
                ++nHitsV;
                break;
            case TPC_VIEW_W:
                ++nHitsW;
                break;
            default:
                break;
        }
    }
    const int nViewsWithHits{(nHitsU > 0 ? 1 : 0) + (nHitsV > 0 ? 1 : 0) + (nHitsW > 0 ? 1 : 0)};
    if (nViewsWithHits < 2)
        return STATUS_CODE_SUCCESS;
 
    const PfoList *pPfoList{nullptr};
    const VertexList *pVertexList{nullptr};
    if (!m_fromList)
    {
        PandoraContentApi::GetCurrentList(*this, pPfoList);
    }
    else
    {
        PandoraContentApi::GetList(*this, m_vertexListName, pVertexList);
    }

    const ParticleFlowObject *pRecoNeutrino{nullptr};
    const MCParticle *pTrueNeutrino{nullptr};
    float energy{-std::numeric_limits<float>::max()};
    for (const MCParticle *const pMC : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMC))
        {
            if (pMC->GetParentList().empty() && pMC->GetEnergy() > energy)
            {
                pTrueNeutrino = pMC;
                energy = pMC->GetEnergy();
            }
        }
    }

    this->Print(pTrueNeutrino, "");

    if (!m_fromList)
    {
        if (pPfoList)
        {
            for (const ParticleFlowObject *pPfo : *pPfoList)
            {
                if (LArPfoHelper::IsNeutrino(pPfo))
                {
                    pRecoNeutrino = pPfo;
                    break;
                }
            }
        }
    }

    const LArMCParticleHelper::FinalStateDescriptor descriptor(*pCaloHitList, *pMCParticleList);

    float trueInelasticity{1.f};
    float trueLeptonEnergy{0.f};
    float trueHadronicEnergy{0.f};
    float trueNeutrinoEnergy{0.f};
    if (pTrueNeutrino)
    {
        trueNeutrinoEnergy = pTrueNeutrino->GetEnergy();
        // Get the fraction of hadroninc energy in the event
        for (const MCParticle *const pChild : pTrueNeutrino->GetDaughterList())
        {
            const int pdg{std::abs(pChild->GetParticleId())};
            if (pdg >= 11 && pdg <= 16)
                trueLeptonEnergy += pChild->GetEnergy();
            else if (pdg != 22)
                trueHadronicEnergy += pChild->GetEnergy();
        }
        if (trueLeptonEnergy > 0.f && trueNeutrinoEnergy > 0.f)
            trueInelasticity = 1.f - trueLeptonEnergy / trueNeutrinoEnergy;
    }

    const bool recoAvailable{pRecoNeutrino || (pVertexList && !pVertexList->empty())};
    if (recoAvailable && pTrueNeutrino)
    {

        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
        const CartesianVector &recoVertex{m_fromList ? pVertexList->front()->GetPosition() : LArPfoHelper::GetVertex(pRecoNeutrino)->GetPosition()};
        if (m_visualise)
        {
            const CartesianVector tu(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ())));
            const CartesianVector tv(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ())));
            const CartesianVector tw(trueVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ())));

            const CartesianVector ru(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoU(recoVertex.GetY(), recoVertex.GetZ())));
            const CartesianVector rv(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoV(recoVertex.GetY(), recoVertex.GetZ())));
            const CartesianVector rw(recoVertex.GetX(), 0.f, static_cast<float>(transform->YZtoW(recoVertex.GetY(), recoVertex.GetZ())));

            const float du{(ru - tu).GetMagnitude()};
            const float dv{(rv - tv).GetMagnitude()};
            const float dw{(rw - tw).GetMagnitude()};

            std::cout << "delta(u, v, w): (" << du << ", " << dv << "," << dw << ")" << std::endl;

            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tu, "U true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tv, "V true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &tw, "W true vertex", BLUE, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &ru, "U reco vertex", RED, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rv, "V reco vertex", RED, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &rw, "W reco vertex", RED, 2));
        }

        if (m_writeFile && LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, "dune_fd_hd"))
        {
            const CartesianVector delta{recoVertex - trueVertex};
            const float dx{delta.GetX()}, dy{delta.GetY()}, dz{delta.GetZ()}, dr{delta.GetMagnitude()};
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
            //std::cout << "E(nu): " << trueNuEnergy << " True Vertex: " << trueVertex << std::endl;
            const int success{1};
            const int nE{static_cast<int>(descriptor.GetNumElectrons())};
            const int nMu{static_cast<int>(descriptor.GetNumMuons())};
            const int nTau{static_cast<int>(descriptor.GetNumTaus())};
            const int nGamma{static_cast<int>(descriptor.GetNumPhotons())};
            const int nP{static_cast<int>(descriptor.GetNumProtons())};
            const int nVisibleP{static_cast<int>(descriptor.GetNumVisibleProtons())};
            const int nN{static_cast<int>(descriptor.GetNumNeutrons())};
            const int nPiC{static_cast<int>(descriptor.GetNumChargedPions())};
            const int nPi0{static_cast<int>(descriptor.GetNumPiZeros())};
            const int nKC{static_cast<int>(descriptor.GetNumChargedKs())};
            const int nK0{static_cast<int>(descriptor.GetNumKZeros())};
            const int nOther{static_cast<int>(descriptor.GetNumOther())};
            const int isCC{static_cast<int>((nE + nMu +  nTau) > 0)};

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueInelasticity", trueInelasticity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nE", nE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMu", nMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTau", nTau));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nGamma));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nP", nP));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nVisibleP", nVisibleP));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nN", nN));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPi0", nPi0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiC", nPiC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nK0", nK0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nKC", nKC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nOther", nOther));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nHitsU", nHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nHitsV", nHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nHitsW", nHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }
    else if (pTrueNeutrino)
    {
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};

        if (m_writeFile && LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, "dune_fd_hd"))
        {
            const int success{0};
            const float dx{-999.f}, dy{-999.f}, dz{-999.f}, dr{-999.f};
            const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
            const int nE{static_cast<int>(descriptor.GetNumElectrons())};
            const int nMu{static_cast<int>(descriptor.GetNumMuons())};
            const int nTau{static_cast<int>(descriptor.GetNumTaus())};
            const int nGamma{static_cast<int>(descriptor.GetNumPhotons())};
            const int nP{static_cast<int>(descriptor.GetNumProtons())};
            const int nN{static_cast<int>(descriptor.GetNumNeutrons())};
            const int nPiC{static_cast<int>(descriptor.GetNumChargedPions())};
            const int nPi0{static_cast<int>(descriptor.GetNumPiZeros())};
            const int nKC{static_cast<int>(descriptor.GetNumChargedKs())};
            const int nK0{static_cast<int>(descriptor.GetNumKZeros())};
            const int nOther{static_cast<int>(descriptor.GetNumOther())};
            const int isCC{static_cast<int>((nE + nMu +  nTau) > 0)};

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueInelasticity", trueInelasticity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nE", nE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMu", nMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTau", nTau));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nGamma));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nP", nP));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nN", nN));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPi0", nPi0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiC", nPiC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nK0", nK0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nKC", nKC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nOther", nOther));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nHitsU", nHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nHitsV", nHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nHitsW", nHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dx", dx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dy", dy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dz", dz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "dr", dr));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }
#endif

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FromList", m_fromList));
    if (m_fromList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
