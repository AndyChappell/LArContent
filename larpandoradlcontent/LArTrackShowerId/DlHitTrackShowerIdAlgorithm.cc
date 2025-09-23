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
#include <cmath>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlHitTrackShowerIdAlgorithm::DlHitTrackShowerIdAlgorithm() :
    m_caloHitListName(""),
    m_trainingMode(false),
    m_vertexRelative(false),
    m_polarCoords(false),
    m_adcPeak(123.f),
    m_maxAdcFactor(10.f),
    m_maxSeqLen(4096)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlHitTrackShowerIdAlgorithm::~DlHitTrackShowerIdAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlHitTrackShowerIdAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
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
    static int event{-1};
    ++event;
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;
    const VertexList *pVertexList(nullptr);
    float vx{0.f};
    float vz{0.f};
    if (m_vertexRelative)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
        if (pVertexList->empty())
        {
            std::cout << "DlHitTrackShowerIdAlgorithm: No vertices found in the vertex list" << std::endl;
            return STATUS_CODE_NOT_FOUND;
        }
        const CartesianVector &vertex{m_vertexRelative ? pVertexList->front()->GetPosition() : CartesianVector(0.f, 0.f, 0.f)};
        vx = vertex.GetX();
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
        switch (pCaloHitList->front()->GetHitType())
        {
            case TPC_VIEW_U:
                vz = transform->YZtoU(vertex.GetY(), vertex.GetZ());
                break;
            case TPC_VIEW_V:
                vz = transform->YZtoV(vertex.GetY(), vertex.GetZ());
                break;
            case TPC_VIEW_W:
                vz = transform->YZtoW(vertex.GetY(), vertex.GetZ());
                break;
            default:
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap, false);
    std::unordered_map<const MCParticle *, const MCParticle *> mcFoldingMap;
    std::unordered_map<const MCParticle *, int> leadingToInstanceMap;
    float xMin{std::numeric_limits<float>::max()}, xMax{std::numeric_limits<float>::lowest()};
    float zMin{std::numeric_limits<float>::max()}, zMax{std::numeric_limits<float>::lowest()};
    float rMax{0};
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const MCParticle *const pLeading{LArMCParticleHelper::GetLeadingEMParticle(pMC)};
        mcFoldingMap[pMC] = pLeading ? pLeading : pMC;
        leadingToInstanceMap[mcFoldingMap.at(pMC)] = 0;
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            xMin = std::min(pos.GetX(), xMin);
            xMax = std::max(pos.GetX(), xMax);
            zMin = std::min(pos.GetZ(), zMin);
            zMax = std::max(pos.GetZ(), zMax);
            if (m_vertexRelative)
            {
                const float x{pos.GetX() - vx};
                const float z{pos.GetZ() - vz};
                const float r{std::sqrt(x * x + z * z)};
                rMax = std::max(r, rMax);
            }
        }
    }
    int inst{1};
    for (const auto &[pLeading, dummy] : leadingToInstanceMap)
    {
        leadingToInstanceMap[pLeading] = inst;
        ++inst;
    }
    const float xRange{xMax - xMin > 0 ? xMax - xMin : 1.f};
    const float zRange{zMax - zMin > 0 ? zMax - zMin : 1.f};

    std::unordered_map<Category, CaloHitList> categoryToHitsMap;
    std::unordered_map<const CaloHit *, std::vector<Category>> hitToCategoriesMap;
    std::unordered_map<const CaloHit *, int> hitToInstanceMap;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC)
        {
            for (const CaloHit *const pCaloHit : caloHitList)
                hitToInstanceMap[pCaloHit] = leadingToInstanceMap.at(mcFoldingMap.at(pMC));
            if (this->IsDiffuse(pLArMC))
            {
                categoryToHitsMap[DIFFUSE].insert(categoryToHitsMap[DIFFUSE].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(DIFFUSE);
            }
            if (this->IsHip(pLArMC))
            {
                categoryToHitsMap[HIP].insert(categoryToHitsMap[HIP].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(HIP);
            }
            if (this->IsMip(pLArMC))
            {
                categoryToHitsMap[MIP].insert(categoryToHitsMap[MIP].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(MIP);
            }
            if (this->IsMichel(pLArMC))
            {
                categoryToHitsMap[LOW_E].insert(categoryToHitsMap[LOW_E].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(LOW_E);
            }
            if (this->IsShower(pLArMC))
            {
                categoryToHitsMap[SHOWER].insert(categoryToHitsMap[SHOWER].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(SHOWER);
            }
            if (this->IsDelta(dynamic_cast<const LArMCParticle *>(mcFoldingMap.at(pLArMC))))
            {
                CaloHitList particleOwnedHits, parentOwnedHits;
                const LArMCParticle *const pParent{dynamic_cast<const LArMCParticle *>(this->AllocateHitOwner(pLArMC, mcFoldingMap, mcToHitsMap,
                    particleOwnedHits, parentOwnedHits))};
                if (!particleOwnedHits.empty())
                {
                    categoryToHitsMap[LOW_E].insert(categoryToHitsMap[LOW_E].end(), particleOwnedHits.begin(), particleOwnedHits.end());
                    for (const CaloHit *const pCaloHit : particleOwnedHits)
                        hitToCategoriesMap[pCaloHit].emplace_back(LOW_E);
                }
                if (!parentOwnedHits.empty())
                {
                    if (this->IsMip(pParent))
                    {
                        categoryToHitsMap[MIP].insert(categoryToHitsMap[MIP].end(), parentOwnedHits.begin(), parentOwnedHits.end());
                        for (const CaloHit *const pCaloHit : parentOwnedHits)
                            hitToCategoriesMap[pCaloHit].emplace_back(MIP);
                    }
                    else
                    {
                        categoryToHitsMap[HIP].insert(categoryToHitsMap[HIP].end(), parentOwnedHits.begin(), parentOwnedHits.end());
                        for (const CaloHit *const pCaloHit : parentOwnedHits)
                            hitToCategoriesMap[pCaloHit].emplace_back(HIP);
                    }
                }
            }
        }
    }
    if (!categoryToHitsMap[DIFFUSE].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[DIFFUSE], "diffuse", CYAN));
    }
    if (!categoryToHitsMap[HIP].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[HIP], "hip", BLACK));
    }
    if (!categoryToHitsMap[MIP].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[MIP], "mip", GRAY));
    }
    if (!categoryToHitsMap[LOW_E].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[LOW_E], "low_e", BLUE));
    }
    if (!categoryToHitsMap[SHOWER].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[SHOWER], "shower", RED));
    }

/*    std::unordered_map<int, CaloHitList> hitInstances;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        for (const CaloHit *const pCaloHit : caloHitList)
            hitInstances[hitToInstanceMap.at(pCaloHit)].emplace_back(pCaloHit);
    }
    for (const auto &[key, hits] : hitInstances)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, "instance_" + std::to_string(key), AUTOITER));
    }*/

    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    categoryToHitsMap.clear();

    // Collect hits and respective categories into feature vectors for training
    FloatVector xx, zz;
    FloatVector rr, cosTheta, sinTheta;
    FloatVector vv;
    if (!m_polarCoords)
    {
        xx.resize(hitToCategoriesMap.size(), 0);
        zz.resize(hitToCategoriesMap.size(), 0);
    }
    else
    {
        rr.resize(hitToCategoriesMap.size(), 0);
        cosTheta.resize(hitToCategoriesMap.size(), 0);
        sinTheta.resize(hitToCategoriesMap.size(), 0);
    }
    vv.resize(hitToCategoriesMap.size(), 0);
    FloatVector adc(hitToCategoriesMap.size(), 0), width(hitToCategoriesMap.size(), 0);
    IntVector semanticLabel(hitToCategoriesMap.size(), UNINITIALISED);
    IntVector instanceLabel(hitToCategoriesMap.size(), UNINITIALISED);

    size_t i{0};
    for (const auto &[pCaloHit, categories] : hitToCategoriesMap)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        const HitType view{pCaloHit->GetHitType()};
        switch (view)
        {
            case TPC_VIEW_U:
                vv[i] = 1 / 3.f;
                break;
            case TPC_VIEW_V:
                vv[i] = 2 / 3.f;
                break;
            case TPC_VIEW_W:
                vv[i] = 3 / 3.f;
                break;
            default:
                return STATUS_CODE_NOT_FOUND;
        }
        if (!m_polarCoords)
        {
            xx[i] = (x - xMin) / xRange;
            zz[i] = (z - zMin) / zRange;
        }
        else
        {
            const float dx{x - vx};
            const float dz{z - vz};
            const float r{std::sqrt(dx * dx + dz * dz)};
            if (r > 0.f)
            {
                rr[i] = r / rMax;
                // Scale to [0,1] to aid training
                cosTheta[i] = 0.5f * (1.f + dx / r);
                sinTheta[i] = 0.5f * (1.f + dz / r);
            }
            else
            {
                rr[i] = 0.f;
                // Scale to [0,1] to aid training
                cosTheta[i] = 0.5f * (1.f + 1.f);
                sinTheta[i] = 0.5f * (1.f + 0.f);
            }
        }
        adc[i] = std::log(1.f + std::min(pCaloHit->GetInputEnergy(), m_adcPeak * m_maxAdcFactor)) / std::log(1.f + m_adcPeak * m_maxAdcFactor);
        width[i] = pCaloHit->GetCellSize1() / xRange;
        for (const auto &category : categories)
        {
            // ATTN: This assumes that the categories are ordered such that MIP < HIP < SHOWER < LOW_E < DIFFUSE
            if (semanticLabel[i] == UNINITIALISED)
            {
                semanticLabel[i] = static_cast<int>(category);
            }
            else
            {
                switch (category)
                {
                    case DIFFUSE:
                        // Allow diffuse to override shower tags
                        if (semanticLabel[i] == SHOWER)
                            semanticLabel[i] = DIFFUSE;
                        break;
                    case HIP:
                    case MIP:
                        // HIP and MIP override everything
                        semanticLabel[i] = category;
                        break;
                    case LOW_E:
                        // Allow low energy to override shower and diffuse tags
                        if (semanticLabel[i] > static_cast<int>(HIP))
                            semanticLabel[i] = LOW_E;
                        break;
                    case SHOWER:
                        // Only tag a hit as shower if it is not initialised
                        break;
                    default:
                        break;
                }
            }
        }
        categoryToHitsMap[Category(semanticLabel[i])].emplace_back(pCaloHit);
        instanceLabel[i] = hitToInstanceMap.at(pCaloHit);
        ++i;
    }

    if (!categoryToHitsMap[DIFFUSE].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[DIFFUSE], "diffuse", CYAN));
    }
    if (!categoryToHitsMap[HIP].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[HIP], "hip", BLACK));
    }
    if (!categoryToHitsMap[MIP].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[MIP], "mip", GRAY));
    }
    if (!categoryToHitsMap[LOW_E].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[LOW_E], "low_e", BLUE));
    }
    if (!categoryToHitsMap[SHOWER].empty())
    {
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &categoryToHitsMap[SHOWER], "shower", RED));
    }

//    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    if (!adc.empty())
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "event", event));
        if (!m_polarCoords)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "view", &vv));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "x", &xx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "z", &zz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "adc", &adc));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "width", &width));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "semantic_label", &semanticLabel));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "instance_label", &instanceLabel));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
        }
        else
        {
            // Sort the vectors by r
            std::vector<size_t> indices(adc.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(), [&rr](size_t a, size_t b) { return rr[a] < rr[b]; });

            auto reorder = [&](auto &vec)
            {
                using VecType = std::decay_t<decltype(vec)>;
                VecType sorted(vec.size());
                for (size_t a = 0; a < indices.size(); ++a)
                    sorted[a] = vec[indices[a]];
                vec.swap(sorted);
            };

            reorder(vv);
            reorder(rr);
            reorder(cosTheta);
            reorder(sinTheta);
            reorder(adc);
            reorder(width);
            reorder(semanticLabel);
            reorder(instanceLabel);

            size_t totalSize{rr.size()};
            for (size_t s = 0; s < totalSize; s += m_maxSeqLen)
            {
                size_t end{std::min(s + m_maxSeqLen, totalSize)};
                FloatVector vvChunk(vv.begin() + s, vv.begin() + end);
                FloatVector rrChunk(rr.begin() + s, rr.begin() + end);
                FloatVector cosThetaChunk(cosTheta.begin() + s, cosTheta.begin() + end);
                FloatVector sinThetaChunk(sinTheta.begin() + s, sinTheta.begin() + end);
                FloatVector adcChunk(adc.begin() + s, adc.begin() + end);
                FloatVector widthChunk(width.begin() + s, width.begin() + end);
                IntVector semanticLabelChunk(semanticLabel.begin() + s, semanticLabel.begin() + end);
                IntVector instanceLabelChunk(instanceLabel.begin() + s, instanceLabel.begin() + end);

                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "view", &vvChunk));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "r", &rrChunk));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "cosTheta", &cosThetaChunk));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "sinTheta", &sinThetaChunk));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "adc", &adcChunk));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "width", &widthChunk));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "semantic_label", &semanticLabelChunk));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "instance_label", &instanceLabelChunk));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
            }
        }
    }

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
    if (parentList.empty() || mcFoldingMap.find(parentList.front()) == mcFoldingMap.end())
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRelative", m_vertexRelative));
    if (m_vertexRelative)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PolarCoords", m_polarCoords));
        if (m_polarCoords)
        {
            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeqLen", m_maxSeqLen));
        }
    }
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AdcPeak", m_adcPeak));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAdcFactor", m_maxAdcFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
