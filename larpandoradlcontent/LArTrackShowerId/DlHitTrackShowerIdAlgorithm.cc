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
    m_trainingMode(false),
    m_adcNormalization(1.f / 123.f)
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
    float xMin{std::numeric_limits<float>::max()}, xMax{std::numeric_limits<float>::lowest()};
    float zMin{std::numeric_limits<float>::max()}, zMax{std::numeric_limits<float>::lowest()};
    float rMax{0};
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const MCParticle *const pLeading{LArMCParticleHelper::GetLeadingEMParticle(pMC)};
        mcFoldingMap[pMC] = pLeading ? pLeading : pMC;
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
    const float xRange{xMax - xMin > 0 ? xMax - xMin : 1.f};
    const float zRange{zMax - zMin > 0 ? zMax - zMin : 1.f};

    std::unordered_map<Category, CaloHitList> categoryToHitsMap;
    std::unordered_map<const CaloHit *, std::vector<Category>> hitToCategoriesMap;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC)
        {
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

//    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    categoryToHitsMap.clear();

    // Collect hits and respective categories into feature vectors for training
    FloatVector xx, zz;
    FloatVector rr, cosTheta, sinTheta;
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
    FloatVector adc(hitToCategoriesMap.size(), 0), width(hitToCategoriesMap.size(), 0);
    IntVector label(hitToCategoriesMap.size(), UNINITIALISED);

    size_t i{0};
    for (const auto &[pCaloHit, categories] : hitToCategoriesMap)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
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
                cosTheta[i] = dx / r;
                sinTheta[i] = dz / r;
            }
            else
            {
                rr[i] = 0.f;
                cosTheta[i] = 1.f;
                sinTheta[i] = 0.f;
            }
        }
        adc[i] = pCaloHit->GetInputEnergy() * m_adcNormalization;
        width[i] = pCaloHit->GetCellSize1() / xRange;
        for (const auto &category : categories)
        {
            // ATTN: This assumes that the categories are ordered such that MIP < HIP < SHOWER < LOW_E < DIFFUSE
            if (label[i] == UNINITIALISED)
            {
                label[i] = static_cast<int>(category);
            }
            else
            {
                switch (category)
                {
                    case DIFFUSE:
                        // Allow diffuse to override shower tags
                        if (label[i] == SHOWER)
                            label[i] = DIFFUSE;
                        break;
                    case HIP:
                    case MIP:
                        // HIP and MIP override everything
                        label[i] = category;
                        break;
                    case LOW_E:
                        // Allow low energy to override shower and diffuse tags
                        if (label[i] > static_cast<int>(HIP))
                            label[i] = LOW_E;
                        break;
                    case SHOWER:
                        // Only tag a hit as shower if it is not initialised
                        break;
                    default:
                        break;
                }
            }
        }
        categoryToHitsMap[Category(label[i])].emplace_back(pCaloHit);
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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "x", &xx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "z", &zz));
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

            reorder(rr);
            reorder(cosTheta);
            reorder(sinTheta);
            reorder(adc);
            reorder(width);
            reorder(label);

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "r", &rr));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "cosTheta", &cosTheta));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "sinTheta", &sinTheta));
        }
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "adc", &adc));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "width", &width));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "label", &label));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
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
    }
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AdcNormalization", m_adcNormalization));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
