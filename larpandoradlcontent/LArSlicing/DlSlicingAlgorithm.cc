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

    if (m_visualize)
        this->VisualizeSlices(sliceToHitsMap);

    this->PopulateRootTree(sliceToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSlicingAlgorithm::Infer()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlSlicingAlgorithm::PopulateRootTree(const LArSliceHelper::SliceHitsMap &mcSlices) const
{
    static int event{-1};
    ++event;
    int sliceId{0};
    const LArTransformationPlugin *const pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    // Fill the training set with slice information
    FloatVector xx, zz;
    IntVector vv, ss, cc;
    float xMin{std::numeric_limits<float>::max()}, xMax{std::numeric_limits<float>::lowest()};
    float uMin{std::numeric_limits<float>::max()}, uMax{std::numeric_limits<float>::lowest()};
    float vMin{std::numeric_limits<float>::max()}, vMax{std::numeric_limits<float>::lowest()};
    float wMin{std::numeric_limits<float>::max()}, wMax{std::numeric_limits<float>::lowest()};

    for (const auto &[pSlice, sliceHits] : mcSlices)
    {
        for (const CaloHit *const pCaloHit : sliceHits)
        {
            xMin = std::min(xMin, pCaloHit->GetPositionVector().GetX());
            xMax = std::max(xMax, pCaloHit->GetPositionVector().GetX());

            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    uMin = std::min(uMin, pCaloHit->GetPositionVector().GetZ());
                    uMax = std::max(uMax, pCaloHit->GetPositionVector().GetZ());
                    break;
                case TPC_VIEW_V:
                    vMin = std::min(vMin, pCaloHit->GetPositionVector().GetZ());
                    vMax = std::max(vMax, pCaloHit->GetPositionVector().GetZ());
                    break;
                case TPC_VIEW_W:
                    wMin = std::min(wMin, pCaloHit->GetPositionVector().GetZ());
                    wMax = std::max(wMax, pCaloHit->GetPositionVector().GetZ());
                    break;
                default:
                    break;
            }
        }
    }
    float xRange{xMax - xMin > 0 ? xMax - xMin : 1.0f};

    for (const auto &[pSlice, sliceHits] : mcSlices)
    {
        const auto &[tpcId, pMC] = pSlice;
        const int pdg{std::abs(pMC->GetParticleId())};
        int view{0};

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
            view = caloHits.front()->GetHitType();
            float zMin{0}, zMax{0};
            switch (view)
            {
                case TPC_VIEW_U:
                    zMin = uMin;
                    zMax = uMax;
                    break;
                case TPC_VIEW_V:
                    zMin = vMin;
                    zMax = vMax;
                    break;
                case TPC_VIEW_W:
                    zMin = wMin;
                    zMax = wMax;
                    break;
                default:
                    break;
            }
            float zRange{zMax - zMin > 0 ? zMax - zMin : 1.0f};

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
            for (const CaloHit *const pCaloHit : caloHits)
            {
                xx.emplace_back((pCaloHit->GetPositionVector().GetX() - xMin) / xRange);
                zz.emplace_back((pCaloHit->GetPositionVector().GetZ() - zMin) / zRange);
                vv.emplace_back(view - TPC_VIEW_U);
                cc.emplace_back(pCaloHit->GetPositionVector() == matchedVertex ? 1 : 0);
                ss.emplace_back(sliceId);
            }
        }

        ++sliceId;
    }

    if (xx.empty())
        return;
    // Don't worry about view separation here, we can handle it easily in HDF5 format conversion if we need to
    // If we do need to split by view, be sure to organise for inference
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "xx", &xx));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "zz", &zz));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "view", &vv));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "cp_label", &cc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "slice_label", &ss));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
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
