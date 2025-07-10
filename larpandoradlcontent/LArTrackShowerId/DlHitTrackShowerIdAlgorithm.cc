/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include <chrono>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlHitTrackShowerIdAlgorithm::DlHitTrackShowerIdAlgorithm() :
    m_caloHitListName(""),
    m_trainingMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlHitTrackShowerIdAlgorithm::~DlHitTrackShowerIdAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::PrepareTrainingSample()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap, false);
    std::unordered_map<const MCParticle *, const MCParticle *> mcFoldingMap;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const MCParticle *const pLeading{LArMCParticleHelper::GetLeadingEMParticle(pMC)};
        mcFoldingMap[pMC] = pLeading ? pLeading : pMC;
    }

    std::map<Category, CaloHitList> categoryToHitsMap;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC)
        {
            if (this->IsDiffuse(pLArMC))
            {
                categoryToHitsMap[DIFFUSE].insert(categoryToHitsMap[DIFFUSE].end(), caloHitList.begin(), caloHitList.end());
            }
            if (this->IsHip(pLArMC))
            {
                categoryToHitsMap[HIP].insert(categoryToHitsMap[HIP].end(), caloHitList.begin(), caloHitList.end());
            }
            if (this->IsMip(pLArMC))
            {
                categoryToHitsMap[MIP].insert(categoryToHitsMap[MIP].end(), caloHitList.begin(), caloHitList.end());
            }
            if (this->IsMichel(pLArMC))
            {
                categoryToHitsMap[LOW_E].insert(categoryToHitsMap[LOW_E].end(), caloHitList.begin(), caloHitList.end());
            }
            if (this->IsShower(pLArMC))
            {
                categoryToHitsMap[SHOWER].insert(categoryToHitsMap[SHOWER].end(), caloHitList.begin(), caloHitList.end());
            }
            if (this->IsDelta(dynamic_cast<const LArMCParticle *>(mcFoldingMap.at(pLArMC))))
            {
                CaloHitList particleOwnedHits, parentOwnedHits;
                const LArMCParticle *const pParent{dynamic_cast<const LArMCParticle *>(this->AllocateHitOwner(pLArMC, mcFoldingMap, mcToHitsMap,
                    particleOwnedHits, parentOwnedHits))};
                if (!particleOwnedHits.empty())
                    categoryToHitsMap[LOW_E].insert(categoryToHitsMap[LOW_E].end(), particleOwnedHits.begin(), particleOwnedHits.end());
                if (!parentOwnedHits.empty())
                {
                    if (this->IsMip(pParent))
                        categoryToHitsMap[MIP].insert(categoryToHitsMap[MIP].end(), parentOwnedHits.begin(), parentOwnedHits.end());
                    else
                        categoryToHitsMap[HIP].insert(categoryToHitsMap[HIP].end(), parentOwnedHits.begin(), parentOwnedHits.end());
                }
            }
        }
    }

    if (!categoryToHitsMap[DIFFUSE].empty())
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[DIFFUSE], "diffuse", CYAN));
    }
    if (!categoryToHitsMap[HIP].empty())
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[HIP], "hip", BLACK));
    }
    if (!categoryToHitsMap[MIP].empty())
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[MIP], "mip", GRAY));
    }
    if (!categoryToHitsMap[LOW_E].empty())
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[LOW_E], "low_e", BLUE));
    }
    if (!categoryToHitsMap[SHOWER].empty())
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[SHOWER], "shower", RED));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::Infer()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsDelta(const LArMCParticle *pLArMC) const
{
    if (!pLArMC)
        return false;
    if (std::abs(pLArMC->GetParticleId()) != E_MINUS)
        return false;

    const MCProcess process{pLArMC->GetProcess()};
    switch (process)
    {
        case MC_PROC_E_IONI:
        case MC_PROC_MU_IONI:
        case MC_PROC_HAD_IONI:
        case MC_PROC_ION_IONI:
            return true;
        default:
            return false;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsDiffuse(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const MCProcess process{pLArMC->GetProcess()};
    if (process == MC_PROC_COMPT)
    {
        if (LArMCParticleHelper::IsDescendentOf(pLArMC, NEUTRON))
            return true;
        if (LArMCParticleHelper::IsDescendentOf(pLArMC, E_MINUS))
            return false;
        else
            return true;
    }
    else
    {
        return false;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsHip(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const int pdg{std::abs(pLArMC->GetParticleId())};
    return (pdg > 1e9 || pdg == PROTON || pdg == K_PLUS || pdg == SIGMA_PLUS || pdg == SIGMA_MINUS || pdg == HYPERON_MINUS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsMip(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const int pdg{std::abs(pLArMC->GetParticleId())};
    return (pdg == MU_MINUS || pdg == PI_PLUS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsMichel(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const MCProcess process{pLArMC->GetProcess()};
    const int pdg{std::abs(pLArMC->GetParticleId())};
    if (process != MC_PROC_DECAY || pdg != E_MINUS)
        return false;
    const MCParticleList &parentList{pLArMC->GetParentList()};
    if (parentList.empty())
        return false;
    const int parentPdg{std::abs(parentList.front()->GetParticleId())};

    // Here we will consider Michel electrons to be those that come from muon decays or pion decays, given the topological similarity
    return (parentPdg == MU_MINUS || parentPdg == PI_PLUS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsShower(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const MCProcess process{pLArMC->GetProcess()};
    const int pdg{std::abs(pLArMC->GetParticleId())};
    return ((pdg == E_MINUS && process != MC_PROC_DECAY) || (pdg == PHOTON && process != MC_PROC_COMPT));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArMCParticle *DlHitTrackShowerIdAlgorithm::AllocateHitOwner(const pandora::MCParticle *const pParticle, const MCFoldingMap &mcFoldingMap,
    const LArMCParticleHelper::MCContributionMap &mcHitMap, CaloHitList &particleOwnedHits, CaloHitList &parentOwnedHits) const
{
    const MCParticleList &parentList{mcFoldingMap.at(pParticle)->GetParentList()};
    const LArMCParticle *pParent{parentList.empty() ? nullptr : dynamic_cast<const LArMCParticle *>(mcFoldingMap.at(parentList.front()))};
    if (pParent == pParticle)
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    while (mcHitMap.find(pParent) == mcHitMap.end())
    {
        pParent = pParent->GetParentList().empty() ? nullptr : dynamic_cast<const LArMCParticle *>(mcFoldingMap.at(pParent->GetParentList().front()));
        if (!pParent)
            return nullptr;
    }
    if (!(this->IsHip(dynamic_cast<const LArMCParticle *>(pParent)) || this->IsMip(dynamic_cast<const LArMCParticle *>(pParent))))
        return nullptr;
    const CaloHitList &particleHits{mcHitMap.at(pParticle)};
    CaloHitList parentHits(mcHitMap.at(pParent).begin(), mcHitMap.at(pParent).end());
    CartesianVector centroid(0, 0, 0);
    LArPcaHelper::EigenValues eigenValues(0, 0, 0);
    LArPcaHelper::EigenVectors eigenVectors;
    LArPcaHelper::RunPca(parentHits, centroid, eigenValues, eigenVectors);
    FloatVector projections;
    for (const CaloHit *const pCaloHit : parentHits)
    {
        const CartesianVector dir{pCaloHit->GetPositionVector() - centroid};
        projections.emplace_back(dir.GetDotProduct(eigenVectors[0]));
    }
    std::vector<size_t> indices(projections.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&projections](size_t i, size_t j) { return projections[i] < projections[j]; });

    CaloHitVector parentHitsVector(parentHits.begin(), parentHits.end());
    std::vector<std::pair<size_t, const CaloHit *>> paired;
    for (size_t i = 0; i < indices.size(); ++i)
        paired.emplace_back(indices[i], parentHitsVector[i]);
    std::sort(paired.begin(), paired.end(), [](const std::pair<size_t, const CaloHit *> &a, const std::pair<size_t, const CaloHit *> &b)
    {
        return a.first > b.first;
    });

    for (const CaloHit *const pCaloHit : particleHits)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        const double xmin{pos.GetX() - 0.5 * pCaloHit->GetCellSize1()}, xmax{pos.GetX() + 0.5 * pCaloHit->GetCellSize1()};
        const double zmin{pos.GetZ() - 0.5 * pCaloHit->GetCellSize0()}, zmax{pos.GetZ() + 0.5 * pCaloHit->GetCellSize0()};
        bool hasInterveningHit{false};
        for (size_t i = 0; i < paired.size() - 1; ++i)
        {
            const CaloHit *const pParentHit1{paired[i].second};
            const CartesianVector &pos1{pParentHit1->GetPositionVector()};
            const CaloHit *const pParentHit2{paired[i + 1].second};
            const CartesianVector &pos2{pParentHit2->GetPositionVector()};
            if (pCaloHit == pParentHit1 || pCaloHit == pParentHit2)
                continue;

            double entry{0.}, exit{1.};
            auto check_axis = [&](double p1, double p2, double minB, double maxB)
            {
                double t1{(minB - p1) / (p2 - p1)};
                double t2{(maxB - p1) / (p2 - p1)};

                if (t1 > t2)
                    std::swap(t1, t2);

                entry = std::max(entry, t1);
                exit = std::min(exit, t2);

                return entry <= exit;
            };

            if (check_axis(pos1.GetX(), pos2.GetX(), xmin, xmax) && check_axis(pos1.GetZ(), pos2.GetZ(), zmin, zmax))
            {
                hasInterveningHit = true;
                break;
            }
        }
        if (!hasInterveningHit)
        {
            particleOwnedHits.emplace_back(pCaloHit);
        }
        else
        {
            parentOwnedHits.emplace_back(pCaloHit);
        }
    }

    return pParent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
