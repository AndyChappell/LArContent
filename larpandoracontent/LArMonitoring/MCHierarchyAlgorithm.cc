/**
 *  @file   larpandoracontent/LArMonitoring/MCHierarchyAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/MCHierarchyAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
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

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    IntVector plane, mcId;
    FloatVector drift, width, channel, adc;
    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        CaloHitList hits3D;
        this->Make3DHits(caloHits, hits3D);
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
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits3D, "3D", AUTOITER));
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

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

void MCHierarchyAlgorithm::Make3DHits(const CaloHitList &hits2D,  CaloHitList &hits3D) const
{
    (void)hits3D;
    if (hits2D.empty())
        return;
    float xMin{std::numeric_limits<float>::max()}, xMax{std::numeric_limits<float>::lowest()};
    for (const CaloHit *const pCaloHit : hits2D)
    {
        const float low{pCaloHit->GetPositionVector().GetX() - 0.5f * pCaloHit->GetCellSize1()};
        const float high{low + pCaloHit->GetCellSize1()};
        if (low < xMin)
            xMin = low;
        if (high > xMax)
            xMax = high;
    }
    if (xMin == xMax)
    {
        xMin -= std::numeric_limits<float>::epsilon();
        xMax += std::numeric_limits<float>::epsilon();
    }
    const int N{static_cast<int>(std::ceil((xMax - xMin) / 0.5f))};

    std::map<int, CaloHitVector> uHits, vHits, wHits;
    for (const CaloHit *const pCaloHit : hits2D)
    {
        const float low{pCaloHit->GetPositionVector().GetX() - 0.5f * pCaloHit->GetCellSize1()};
        const float high{low + pCaloHit->GetCellSize1()};
        const int bin0{static_cast<int>(std::floor((low - xMin) / 0.5f))};
        const int bin1{static_cast<int>(std::floor((high - xMin) / 0.5f))};
        for (int b = bin0; b <= bin1; ++b)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    uHits[b].emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    vHits[b].emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    wHits[b].emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
    }

    std::map<int, std::map<int, std::map<int, float>>> matrix;
    CaloHitSet usedHitsU, usedHitsV, usedHitsW;
    // In each bin, identify the best combination of hits based on chi2, reject chi2 > 1.5 (maybe 6 due to dof)
    // Need to track used hits, their best chi2, the 3D pos and the hits that accompany them
    for (int b = 0; b < N; ++b)
    {
        if (uHits[b].empty() || vHits[b].empty() || wHits[b].empty())
            continue;
        matrix.clear();

        // Index using [] because it's likely we'll have gaps in the binning which we can just skip
        for (size_t i = 0; i < uHits[b].size(); ++i)
        {
            const CaloHit *const pCaloHitU{uHits[b].at(i)};
            const CartesianVector &posU{pCaloHitU->GetPositionVector()};
            for (size_t j = 0; j < vHits[b].size(); ++j)
            {
                const CaloHit *const pCaloHitV{vHits[b].at(j)};
                const CartesianVector &posV{pCaloHitV->GetPositionVector()};
                for (size_t k = 0; k < wHits[b].size(); ++k)
                {
                    const CaloHit *const pCaloHitW{wHits[b].at(k)};
                    const CartesianVector &posW{pCaloHitW->GetPositionVector()};
                    float chi2{std::numeric_limits<float>::max()};
                    CartesianVector pos3D(0, 0, 0);
                    LArGeometryHelper::MergeThreeWidePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, posU, posV, posW,
                        0.5f * pCaloHitU->GetCellSize1(), 0.5f * pCaloHitV->GetCellSize1(), 0.5f * pCaloHitW->GetCellSize1(), pos3D, chi2);
                    matrix[i][j][k] = chi2;
                }
            }
        }

        // Check what's left over, maybe make 2 hit 3D hits from the leftover cases, perhaps checking proximity
        // to the existing confirmed hits as a way to ensure some kind of consistency (enfore closest approach of e.g. < 5 cm)
        std::vector<CandidateHit> candidateHits;
        for (size_t i = 0; i < uHits[b].size(); ++i)
        {
            for (size_t j = 0; j < vHits[b].size(); ++j)
            {
                for (size_t k = 0; k < wHits[b].size(); ++k)
                    candidateHits.emplace_back(std::make_tuple(i, j, k, matrix[i][j][k]));
            }
        }
        std::sort(candidateHits.begin(), candidateHits.end(), [] (const CandidateHit &hit1, const CandidateHit &hit2)
            { return std::get<3>(hit1) < std::get<3>(hit2); });
        // Note these should be promoted to broader scope, because in principle hits can turn up in multiple bins
        std::vector<CandidateHit> confirmedHits;
        for (const CandidateHit &candidateHit : candidateHits)
        {
            if (usedHitsU.find(uHits[b].at(std::get<0>(candidateHit))) != usedHitsU.end())
                continue;
            if (usedHitsV.find(vHits[b].at(std::get<1>(candidateHit))) != usedHitsV.end())
                continue;
            if (usedHitsW.find(wHits[b].at(std::get<2>(candidateHit))) != usedHitsW.end())
                continue;
            if (std::get<3>(candidateHit) > 6)
                continue;
            usedHitsU.insert(uHits[b].at(std::get<0>(candidateHit)));
            usedHitsV.insert(vHits[b].at(std::get<1>(candidateHit)));
            usedHitsW.insert(wHits[b].at(std::get<2>(candidateHit)));
            confirmedHits.emplace_back(candidateHit);
        }

        for (const CandidateHit &hit : confirmedHits)
        {
            const CaloHit *pCaloHit3D{nullptr};
            const CaloHit *uHit{uHits[b].at(std::get<0>(hit))};
            const CaloHit *vHit{vHits[b].at(std::get<1>(hit))};
            const CaloHit *wHit{wHits[b].at(std::get<2>(hit))};
            float chi2{0.f};
            CartesianVector pos3D(0, 0, 0);
            LArGeometryHelper::MergeThreeWidePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W,
                uHit->GetPositionVector(), vHit->GetPositionVector(), wHit->GetPositionVector(), 0.5f * uHit->GetCellSize1(),
                0.5f * vHit->GetCellSize1(), 0.5f * wHit->GetCellSize1(), pos3D, chi2);

            PandoraContentApi::CaloHit::Parameters parameters;
            parameters.m_hitType = TPC_3D;
            const CaloHit *const pCaloHit2D(wHit);
            parameters.m_positionVector = pos3D;
            parameters.m_pParentAddress = static_cast<const void *>(pCaloHit2D);
            parameters.m_cellThickness = pCaloHit2D->GetCellThickness();
            parameters.m_cellGeometry = RECTANGULAR;
            parameters.m_cellSize0 = pCaloHit2D->GetCellLengthScale();
            parameters.m_cellSize1 = pCaloHit2D->GetCellLengthScale();
            parameters.m_cellNormalVector = pCaloHit2D->GetCellNormalVector();
            parameters.m_expectedDirection = pCaloHit2D->GetExpectedDirection();
            parameters.m_nCellRadiationLengths = pCaloHit2D->GetNCellRadiationLengths();
            parameters.m_nCellInteractionLengths = pCaloHit2D->GetNCellInteractionLengths();
            parameters.m_time = pCaloHit2D->GetTime();
            parameters.m_inputEnergy = pCaloHit2D->GetInputEnergy();
            parameters.m_mipEquivalentEnergy = pCaloHit2D->GetMipEquivalentEnergy();
            parameters.m_electromagneticEnergy = pCaloHit2D->GetElectromagneticEnergy();
            parameters.m_hadronicEnergy = pCaloHit2D->GetHadronicEnergy();
            parameters.m_isDigital = pCaloHit2D->IsDigital();
            parameters.m_hitRegion = pCaloHit2D->GetHitRegion();
            parameters.m_layer = pCaloHit2D->GetLayer();
            parameters.m_isInOuterSamplingLayer = pCaloHit2D->IsInOuterSamplingLayer();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));
            hits3D.emplace_back(pCaloHit3D);
        }
    }
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