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
    m_hitsTreeName{"hits"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCHierarchyAlgorithm::Run()
{
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
        {
            std::cout << "Found root " << pMC->GetParticleId() << std::endl;
            pRoot = pMC;
        }
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
    std::cout << "Map Length: " << m_mcToHitsMap.size() << " List Length: " << mcList.size() << std::endl;

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
        while (!pParent->GetParentList().empty())
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
        // Want event level tree for things like CC/NC, flavour, num final state particles and possibly final state track/shower count
        // Need a separate tree for the hits, x, z, adc, and also need to consider creating 3D particles
        // with x, y, z, adc from the MC
        const int mcId{mcToIdMap.at(pMC)};
        const CartesianVector &mom{pMC->GetMomentum()};
        FloatVector momVec{mom.GetX(), mom.GetY(), mom.GetZ()};
        const CartesianVector &vtx{pMC->GetParticleId() != 22 ? pMC->GetVertex() : pMC->GetEndpoint()};
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

        // Need to construct arbitrary id for MC so we can have parent and child indices
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "event_id", event));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "mc_id", mcId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "pdg", pMC->GetParticleId()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "energy", pMC->GetEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "mom_vec", &momVec));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "parent_id", parentId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "child_id_vec", &childIdVector));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_total", nHitsU + nHitsV + nHitsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_u", nHitsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_v", nHitsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "n_hits_w", nHitsW));
        // Want a particle tier related to the visible particles
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_vec", &vtxVec));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_x", xVtx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_u", uVtx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_v", vVtx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "vtx_w", wVtx));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_mcTreeName));
    }


/*
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "isCC", isCC));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "isNuE", isNuE));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "isNuMu", isNuMu));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_mcTreeName, "nFinalStateParticles", nFinalStateParticles));
*/

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

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
