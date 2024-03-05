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

    MCParticleList primariesList(primaries.begin(), primaries.end());
    const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primariesList)};

    float trueInelasticity{1.f};
    if (pTrueNeutrino)
    {
        // Get the fraction of hadroninc energy in the event
        float trueLeptonEnergy{0.f};
        for (const MCParticle *const pChild : pTrueNeutrino->GetDaughterList())
        {
            const int pdg{std::abs(pChild->GetParticleId())};
            if (pdg >= 11 && pdg <= 16)
                trueLeptonEnergy += pChild->GetEnergy();
        }
        if (trueLeptonEnergy > 0.f && pTrueNeutrino->GetEnergy() > 0.f)
            trueInelasticity = 1.f - trueLeptonEnergy / pTrueNeutrino->GetEnergy();
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
            const int success{1};
            const int isCC{descriptor.IsCC()};
            const int isQE{descriptor.IsQE()};
            const int isRes{descriptor.IsResonant()};
            const int isDIS{descriptor.IsDIS()};
            const int isCoh{descriptor.IsCoherent()};
            const int isOther{!(isCC || isQE || isRes || isDIS)};
            const int isMu{descriptor.IsMuonNeutrino()};
            const int isElectron{descriptor.IsElectronNeutrino()};
            const int nPiZero{static_cast<int>(descriptor.GetNumPiZero())};
            const int nPiPlus{static_cast<int>(descriptor.GetNumPiPlus())};
            const int nPiMinus{static_cast<int>(descriptor.GetNumPiMinus())};
            const int nPhotons{static_cast<int>(descriptor.GetNumPhotons())};
            const int nProtons{static_cast<int>(descriptor.GetNumProtons())};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueInelasticity", trueInelasticity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isQE", isQE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isRes", isRes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isDIS", isDIS));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCoh", isCoh));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isOther", isOther));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMu", isMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isElectron", isElectron));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiZero", nPiZero));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiPlus", nPiPlus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiMinus", nPiMinus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nPhotons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nProtons", nProtons));
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
            const int isCC{descriptor.IsCC()};
            const int isQE{descriptor.IsQE()};
            const int isRes{descriptor.IsResonant()};
            const int isDIS{descriptor.IsDIS()};
            const int isCoh{descriptor.IsCoherent()};
            const int isOther{!(isCC || isQE || isRes || isDIS)};
            const int isMu{descriptor.IsMuonNeutrino()};
            const int isElectron{descriptor.IsElectronNeutrino()};
            const int nPiZero{static_cast<int>(descriptor.GetNumPiZero())};
            const int nPiPlus{static_cast<int>(descriptor.GetNumPiPlus())};
            const int nPiMinus{static_cast<int>(descriptor.GetNumPiMinus())};
            const int nPhotons{static_cast<int>(descriptor.GetNumPhotons())};
            const int nProtons{static_cast<int>(descriptor.GetNumProtons())};
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "trueInelasticity", trueInelasticity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isQE", isQE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isRes", isRes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isDIS", isDIS));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCoh", isCoh));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isOther", isOther));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMu", isMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isElectron", isElectron));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiZero", nPiZero));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiPlus", nPiPlus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiMinus", nPiMinus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nPhotons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nProtons", nProtons));
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
