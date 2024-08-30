/**
 *  @file   larpandoracontent/LArMonitoring/MCHierarchyAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/MCHierarchyAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

MCHierarchyAlgorithm::MCHierarchyAlgorithm() :
    m_caloHitListName{"CaloHitList2D"},
    m_mcParticleListName{"Input"},
    m_mcFileName{"mc.root"},
    m_mcTreeName{"mc"},
    m_eventFileName{"events.root"},
    m_eventTreeName{"events"},
    m_hitsFileName{"hits.root"},
    m_hitsTreeName{"hits"},
    m_correctionX{0.f},
    m_correctionY{0.f},
    m_correctionZ{0.f},
    m_visualize{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCHierarchyAlgorithm::~MCHierarchyAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_mcTreeName, m_mcFileName, "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_hitsTreeName, m_hitsFileName, "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_eventTreeName, m_eventFileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCHierarchyAlgorithm::Run()
{
    m_mcToHitsMap.clear();
    const CartesianVector geoCorrection(m_correctionX, m_correctionY, m_correctionZ);
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }

    static int event{-1};
    ++event;
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    if (!pMCParticleList || pMCParticleList->empty() || !pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    typedef std::map<const MCParticle *, const MCParticle *> MCToMCMap;
    const MCParticle *pRoot{nullptr};
    MCToMCMap mcToLeadingMap;
    for (const MCParticle *const pMC : *pMCParticleList)
    {
        if ((LArMCParticleHelper::IsNeutrino(pMC) || LArMCParticleHelper::IsTriggeredBeamParticle(pMC)) && pMC->GetParentList().empty())
            pRoot = pMC;
        const MCParticle *pParent{pMC}, *pLeadingEM{nullptr};
        while (!pParent->GetParentList().empty())
        {
            pParent = pParent->GetParentList().front();
            const int pdg{std::abs(pParent->GetParticleId())};
            if (pdg == E_MINUS || pdg == PHOTON)
                pLeadingEM = pParent;
        }
        if (pLeadingEM)
            mcToLeadingMap[pMC] = pLeadingEM;
        else
            mcToLeadingMap[pMC] = pMC;
    }

    if (!pRoot)
        return STATUS_CODE_SUCCESS;
    const CartesianVector trueVtx(pRoot->GetVertex().GetX() + geoCorrection.GetX(), 0, pRoot->GetVertex().GetZ() + geoCorrection.GetZ());
    // Veto events interacting in the beamline instrumentation
    if (trueVtx.GetZ() < -2.f)
        return STATUS_CODE_SUCCESS;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            if (mcToLeadingMap.find(pMCParticle) == mcToLeadingMap.end())
                continue;
            const MCParticle *const pLeading{mcToLeadingMap.at(pMCParticle)};
            m_mcToHitsMap[pLeading].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &exception)
        {
        }
    }
    // The root neutrino will not be added by the preceding code, so add it
    if (m_mcToHitsMap.find(pRoot) == m_mcToHitsMap.end())
        m_mcToHitsMap[pRoot] = CaloHitList();

    MCParticleList mcList;
    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        if (pMC != pRoot)
            mcList.emplace_back(pMC);
        else
            mcList.emplace_front(pMC);
    }

    int i{0};
    typedef std::map<const MCParticle *, int> MCToIndexMap;
    MCToIndexMap mcToIdMap;
    for (const MCParticle *const pMC : mcList)
    {
        mcToIdMap[pMC] = i;
        ++i;
    }
 
    // Determine map between visible MC particles
    MCToMCMap mcParentMap;
    typedef std::map<const MCParticle *, MCParticleList> MCToMCListMap;
    MCToMCListMap mcChildMap;
    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        const MCParticle *pParent{pMC};
        while (pParent->GetParentList().size() == 1)
        {
            pParent = pParent->GetParentList().front();
            if (m_mcToHitsMap.find(pParent) != m_mcToHitsMap.end())
            {
                mcParentMap[pMC] = pParent;
                mcChildMap[pParent].emplace_back(pMC);
                break;
            }
        }
        if (mcParentMap.find(pMC) == mcParentMap.end())
            mcParentMap[pMC] = nullptr;
    }

    // Ensure any particles without children still have an entry in the MC->Child map
    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        if (mcChildMap.find(pMC) == mcChildMap.end())
            mcChildMap[pMC] = MCParticleList();
    }

    const LArTransformationPlugin *const transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        const int mcId{mcToIdMap.at(pMC)};
        const int isTriggeredBeam{LArMCParticleHelper::IsTriggeredBeamParticle(pMC) ? 1 : 0};
        const int isBeamInduced{!isTriggeredBeam && LArMCParticleHelper::IsBeamParticle(pMC) ? 1 : 0};
        const CartesianVector &mom{pMC->GetMomentum()};
        FloatVector momVec{mom.GetX(), mom.GetY(), mom.GetZ()};
        const CartesianVector &vtx{pMC->GetParticleId() != 22 ? pMC->GetVertex() + geoCorrection : pMC->GetEndpoint() + geoCorrection};
        FloatVector vtxVec{vtx.GetX(), vtx.GetY(), vtx.GetZ()};
        const float xVtx{vtx.GetX()};
        const float uVtx{static_cast<float>(transform->YZtoU(vtx.GetY(), vtx.GetZ()))};
        const float vVtx{static_cast<float>(transform->YZtoV(vtx.GetY(), vtx.GetZ()))};
        const float wVtx{static_cast<float>(transform->YZtoW(vtx.GetY(), vtx.GetZ()))};
        const MCParticle *const parent{mcParentMap.at(pMC)};
        const int parentId{parent ? mcToIdMap.at(parent) : -1};
        IntVector childIdVector;
        for (const MCParticle *const pChild : mcChildMap.at(pMC))
            childIdVector.emplace_back(mcToIdMap.at(pChild));

        CaloHitList uHits, vHits, wHits;
        for (const CaloHit *const pCaloHit : caloHits)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    uHits.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    vHits.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    wHits.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
        int nHitsU{static_cast<int>(uHits.size())}, nHitsV{static_cast<int>(vHits.size())}, nHitsW{static_cast<int>(wHits.size())};

        if (m_visualize)
        {
            if (LArMCParticleHelper::IsTriggeredBeamParticle(pMC))
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "trig", RED));
            }
            else if (LArMCParticleHelper::IsBeamParticle(pMC))
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "beam", BLACK));
            }
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "event_id", event));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "mc_id", mcId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "is_triggered_beam", isTriggeredBeam));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "is_beam_induced", isBeamInduced));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "pdg", pMC->GetParticleId()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "energy", pMC->GetEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "mom_vec", &momVec));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "parent_id", parentId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "child_id_vec", &childIdVector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_total", nHitsU + nHitsV + nHitsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_u", nHitsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_v", nHitsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_w", nHitsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_vec", &vtxVec));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_x", xVtx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_u", uVtx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_v", vVtx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_w", wVtx));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_mcTreeName));
    }
    if (m_visualize)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueVtx, "true", BLUE, 3));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    // Consider creating 3D particles with x, y, z, adc from the MC
    IntVector plane, mcId;
    FloatVector drift, width, channel, adc;
    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        for (const CaloHit *const pCaloHit : caloHits)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_W:
                    plane.emplace_back(0);
                    break;
                case TPC_VIEW_V:
                    plane.emplace_back(1);
                    break;
                case TPC_VIEW_U:
                    plane.emplace_back(2);
                    break;
                default:
                    continue;
            }
            mcId.emplace_back(mcToIdMap.at(pMC));
            drift.emplace_back(pCaloHit->GetPositionVector().GetX());
            width.emplace_back(pCaloHit->GetCellSize1());
            channel.emplace_back(pCaloHit->GetPositionVector().GetZ());
            adc.emplace_back(pCaloHit->GetInputEnergy());
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_hitsTreeName, "event_id", event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_hitsTreeName, "plane", &plane));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_hitsTreeName, "mc_id", &mcId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_hitsTreeName, "drift", &drift));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_hitsTreeName, "width", &width));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_hitsTreeName, "channel", &channel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_hitsTreeName, "adc", &adc));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_hitsTreeName));

    const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pRoot)};
    const int isCC{pLArMC->IsCC()};
    const int isNue{std::abs(pLArMC->GetParticleId()) == 12};
    const int isNumu{std::abs(pLArMC->GetParticleId()) == 14};
    const int isNutau{std::abs(pLArMC->GetParticleId()) == 16};
    const int isTestbeam(LArMCParticleHelper::IsTriggeredBeamParticle(pRoot) ? 1 : 0);
    int decaysHadronically{0};
    if (isCC && isNutau)
    {
        for (const MCParticle *const pChild : pLArMC->GetDaughterList())
        {
            const int pdg{std::abs(pChild->GetParticleId())};
            if (pdg == E_MINUS || pdg == MU_MINUS)
            {
                decaysHadronically = 1;
                break;
            }
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "event_id", event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "is_cc", isCC));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "is_nu_e", isNue));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "is_nu_mu", isNumu));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "is_nu_tau", isNutau));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "tau_decays_hadronically", decaysHadronically));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "is_testbeam", isTestbeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName, "n_visible_final_state_particles",
        static_cast<int>(mcChildMap.at(pRoot).size())));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_eventTreeName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCFileName", m_mcFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCTreeName", m_mcTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventFileName", m_eventFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventTreeName", m_eventTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitsFileName", m_hitsFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitsTreeName", m_hitsTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CorrectionX", m_correctionX));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CorrectionY", m_correctionY));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CorrectionZ", m_correctionZ));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
