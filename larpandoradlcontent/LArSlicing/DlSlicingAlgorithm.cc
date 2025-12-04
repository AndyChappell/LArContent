/**
 *  @file   larpandoracontent/LArMonitoring/DlSlicingAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include <list>
#include <unordered_set>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArSliceHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoradlcontent/LArSlicing/DlSlicingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DlSlicingAlgorithm::DlSlicingAlgorithm() :
    m_trainingMode(false),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlSlicingAlgorithm::~DlSlicingAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
    }
    catch (StatusCodeException e)
    {
        std::cout << "DlSlicingAlgorithm: Unable to write to ROOT tree" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::Run()
{
    if (m_trainingMode)
    {
        return this->PrepareTrainingSample();
    }
    else
    {
        return this->Infer();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::PrepareTrainingSample()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap, false);

    LArMCParticleHelper::MCLeadingMap mcToLeadingMap;
    LArMCParticleHelper::GetMCToLeadingMap(mcToHitsMap, mcToLeadingMap);

    LArSliceHelper::SliceHitsMap sliceToHitsMap;
    LArSliceHelper::GetSliceToHitsMap(mcToHitsMap, mcToLeadingMap, sliceToHitsMap);

    CaloHitList backgroundHits;
    LArSliceHelper::FilterSlices({NEUTRON, PHOTON}, sliceToHitsMap, backgroundHits);

    if (m_visualize)
        this->VisualizeSlices(sliceToHitsMap, backgroundHits);

    this->PopulateRootTree(sliceToHitsMap, backgroundHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::Infer()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlSlicingAlgorithm::PopulateRootTree(const LArSliceHelper::SliceHitsMap &mcSlices, const CaloHitList &backgroundHits) const
{
    static int event{-1};
    ++event;
    int sliceId{0};
    const LArTransformationPlugin *const pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    // Fill the training set with per slice information
    for (const auto &[pSlice, sliceHits] : mcSlices)
    {
        const auto &[tpcId, pMC] = pSlice;
        const int pdg{std::abs(pMC->GetParticleId())};
        const int isNeutrino{(pdg == NU_E || pdg == NU_MU || pdg == NU_TAU) ? 1 : 0};
        int view{0};

        FloatVector xx, zz;

        // Separate hits by view
        CaloHitList caloHitsU, caloHitsV, caloHitsW;
        for (const CaloHit *const pCaloHit : sliceHits)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    caloHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    caloHitsV.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    caloHitsW.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }

        for (const CaloHitList &caloHits : {caloHitsU, caloHitsV, caloHitsW})
        {
            if (caloHits.empty())
                continue;
            xx.clear(); zz.clear();
            view = caloHits.front()->GetHitType();
            for (const CaloHit *const pCaloHit : caloHits)
            {
                xx.emplace_back(pCaloHit->GetPositionVector().GetX());
                zz.emplace_back(pCaloHit->GetPositionVector().GetZ());
            }
            CartesianVector matchedVertex(0, 0, 0);
            if (pdg == MU_MINUS)
            {
                CartesianVector trueVertex(0, 0, 0);
                LArVertexHelper::GetProjectedTrueVertex(pTransform, pMC, caloHits.front()->GetHitType(), trueVertex);

                CartesianVector trueDirection{pMC->GetMomentum().GetUnitVector()};
                // Our hit lists here are for the whole slice, but for cosmics, we only want to consider hits from the cosmic muon itself to avoid
                // attaching a vertex to a child particle's hits
                CaloHitList filteredHits;
                this->FilterSliceHitsToCosmic(caloHits, pMC, filteredHits);
                LArVertexHelper::MatchHitToCosmicVertex(pTransform, filteredHits, trueDirection, matchedVertex);
            }
            else
            {
                CartesianVector trueVertex(0, 0, 0);
                LArVertexHelper::GetProjectedTrueVertex(pTransform, pMC, caloHits.front()->GetHitType(), trueVertex);
                LArVertexHelper::MatchHitToVertex(caloHits, trueVertex, matchedVertex);
            }

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "event", event));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "slice_id", sliceId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "is_true_neutrino", isNeutrino));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "is_background", 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "view", view));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "vertex_x", matchedVertex.GetX()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "vertex_z", matchedVertex.GetZ()));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "xx", &xx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "zz", &zz));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
        }

        ++sliceId;
    }

    // Fill the training set with background hit information
    if (!backgroundHits.empty())
    {
        FloatVector xx, zz;

        // Separate hits by view
        CaloHitList caloHitsU, caloHitsV, caloHitsW;
        for (const CaloHit *const pCaloHit : backgroundHits)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    caloHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    caloHitsV.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    caloHitsW.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }

        for (const CaloHitList &caloHits : {caloHitsU, caloHitsV, caloHitsW})
        {
            if (caloHits.empty())
                continue;
            xx.clear(); zz.clear();
            const int view{caloHits.front()->GetHitType()};
            for (const CaloHit *const pCaloHit : caloHits)
            {
                xx.emplace_back(pCaloHit->GetPositionVector().GetX());
                zz.emplace_back(pCaloHit->GetPositionVector().GetZ());
            }
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "event", event));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "slice_id", sliceId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "is_true_neutrino", 0));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "is_background", 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "view", view));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "xx", &xx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "zz", &zz));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlSlicingAlgorithm::FilterSliceHitsToCosmic(const CaloHitList &sliceHits, const MCParticle *const pCosmicMC,
    CaloHitList &cosmicHits) const
{
    for (const CaloHit *const pCaloHit : sliceHits)
    {
        try
        {
            const MCParticle *const pMainMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            if (pMainMC == pCosmicMC)
                cosmicHits.emplace_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
            continue;
        }
    }

    // In the unlikely event that the cosmic ray has no direct hits, return the original slice hits
    if (cosmicHits.empty())
        cosmicHits = sliceHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlSlicingAlgorithm::VisualizeSlices(const LArSliceHelper::SliceHitsMap &mcSlices) const
{
    const CaloHitList dummyList;
    this->VisualizeSlices(mcSlices, dummyList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlSlicingAlgorithm::VisualizeSlices(const LArSliceHelper::SliceHitsMap &mcSlices, const CaloHitList &backgroundHits) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));

    for (const auto &[slice, caloHitList] : mcSlices)
    {
        const MCParticle *const pMC{slice.second};
        CaloHitList caloHitsU, caloHitsV, caloHitsW;
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    caloHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    caloHitsV.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    caloHitsW.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
        const int pdg{std::abs(pMC->GetParticleId())};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsU, "MC(U): " + std::to_string(pdg), RED));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsV, "MC(V): " + std::to_string(pdg), GREEN));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsW, "MC(W): " + std::to_string(pdg), BLUE));
    }

    if (!backgroundHits.empty())
    {
        CaloHitList caloHitsU, caloHitsV, caloHitsW;
        for (const CaloHit *const pCaloHit : backgroundHits)
        {
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    caloHitsU.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    caloHitsV.emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    caloHitsW.emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsU, "Background U", GRAY));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsV, "Background V", GRAY));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitsW, "Background W", GRAY));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
