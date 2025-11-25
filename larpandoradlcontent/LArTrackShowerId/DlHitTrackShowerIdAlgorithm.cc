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
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <unistd.h>

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
    float vx{0.f};
    float vz{0.f};
    if (m_vertexRelative)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, != , this->GetVertexPlanarCoordinates(pCaloHitList->front()->GetHitType(), vx, vz));
    }

    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(*pCaloHitList, mcToHitsMap, false);
    MCFoldingMap mcFoldingMap;
    Bounds bounds;
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        const MCParticle *const pLeading{LArMCParticleHelper::GetLeadingEMParticle(pMC)};
        mcFoldingMap[pMC] = pLeading ? pLeading : pMC;
        this->GetCoordinateExtrema(caloHitList, vx, vz, bounds);
    }

    lar_content::LArMCParticleHelper::MCContributionMap leadingToHitsMap;
    this->FoldToLeading(mcToHitsMap, mcFoldingMap, leadingToHitsMap);
    lar_content::LArMCParticleHelper::MCContributionMap instanceToHitsMap;
    this->ConsolidateInstances(leadingToHitsMap, instanceToHitsMap);

    int inst{1};
    std::unordered_map<const MCParticle *, int> mcToInstanceMap;
    for (const auto &[pMC, dummy] : instanceToHitsMap)
    {
        mcToInstanceMap[pMC] = inst;
        ++inst;
    }

    std::unordered_map<Category, CaloHitList> categoryToHitsMap;
    std::unordered_map<const CaloHit *, std::vector<Category>> hitToCategoriesMap;
    std::unordered_map<const CaloHit *, int> hitToInstanceMap;
    for (const auto &[pMC, caloHitList] : instanceToHitsMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC)
        {
            for (const CaloHit *const pCaloHit : caloHitList)
                hitToInstanceMap[pCaloHit] = mcToInstanceMap.at(pMC);
            if (this->IsHip(pLArMC))
            {
                categoryToHitsMap[HIP].insert(categoryToHitsMap[HIP].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(HIP);
            }
            else if (this->IsMip(pLArMC))
            {
                categoryToHitsMap[MIP].insert(categoryToHitsMap[MIP].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(MIP);
            }
            else if (this->IsDelta(pLArMC))
            {
                categoryToHitsMap[LOW_E].insert(categoryToHitsMap[LOW_E].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(LOW_E);
            }
            else if (this->IsShower(pLArMC))
            {
                categoryToHitsMap[SHOWER].insert(categoryToHitsMap[SHOWER].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(SHOWER);
            }
            else if (this->IsMichel(pLArMC))
            {
                categoryToHitsMap[LOW_E].insert(categoryToHitsMap[LOW_E].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(LOW_E);
            }
            else if (this->IsDiffuse(pLArMC))
            {
                // Don't believe we'll ever get here in the folded scenario, but class diffuse as shower-like just in case
                categoryToHitsMap[SHOWER].insert(categoryToHitsMap[SHOWER].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(SHOWER);
            }
            else
            {
                std::cout << "Warning: Found unclassifiable particle, default to MIP" << std::endl;
                std::cout << "  " << pLArMC->GetParticleId() << " " << pLArMC->GetProcess() << ": " << caloHitList.size() << " hits" << std::endl;
                std::cout << "  Parent PDG code: " << (pLArMC->GetParentList().empty() ? 0 : pLArMC->GetParentList().front()->GetParticleId()) << std::endl;
                categoryToHitsMap[MIP].insert(categoryToHitsMap[MIP].end(), caloHitList.begin(), caloHitList.end());
                for (const CaloHit *const pCaloHit : caloHitList)
                    hitToCategoriesMap[pCaloHit].emplace_back(MIP);
            }
        }
    }

    for (const auto &[key, hits] : categoryToHitsMap)
    {
        std::string categoryName;
        switch (key)
        {
            case MIP:
                categoryName = "MIP";
                break;
            case HIP:
                categoryName = "HIP";
                break;
            case SHOWER:
                categoryName = "SHOWER";
                break;
            case LOW_E:
                categoryName = "LOW_E";
                break;
            default:
                categoryName = "UNINITIALISED";
                break;
        }
        std::cout << categoryName << ": " << hits.size() << " hits" << std::endl;
    }

    /*
    for (const auto &[cat, caloHitList] : categoryToHitsMap)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, std::to_string(cat), AUTOITER));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    for (const auto &[pMC, caloHitList] : instanceToHitsMap)
    {
        const MCParticleList &parentList{pMC->GetParentList()};
        const std::string str{parentList.empty() ? "" : std::to_string(std::abs(parentList.front()->GetParticleId())) + " -> "};
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, str + std::to_string(pMC->GetParticleId()), AUTOITER));
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    */
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
            xx[i] = (x - bounds.xMin) / bounds.xRange;
            zz[i] = (z - bounds.zMin) / bounds.zRange;
        }
        else
        {
            const float dx{x - vx};
            const float dz{z - vz};
            const float r{std::sqrt(dx * dx + dz * dz)};
            if (r > 0.f)
            {
                rr[i] = r / bounds.rMax;
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
        width[i] = pCaloHit->GetCellSize1() / bounds.xRange;
        for (const auto &category : categories)
            semanticLabel[i] = static_cast<int>(category);
        instanceLabel[i] = hitToInstanceMap.at(pCaloHit);
        ++i;
    }

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
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    // Get projected vertex coordinates if needed
    float vx{0.f};
    float vz{0.f};
    if (m_vertexRelative)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, != , this->GetVertexPlanarCoordinates(pCaloHitList->front()->GetHitType(), vx, vz));
    }

    // Get coordinate minima and maxima
    Bounds bounds;
    this->GetCoordinateExtrema(*pCaloHitList, vx, vz, bounds);

    // Prepare input feature vectors
    CaloHitVector caloHitVector;
    FloatVector xx, zz , rr, cosTheta, sinTheta, vv, adc, width;
    this->PopulateInputVectors(*pCaloHitList, bounds, vx, vz, xx, zz, rr, cosTheta, sinTheta, vv, adc, width, caloHitVector);

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));
    // Ensure each run conforms to the maximum sequence length
    size_t totalSize{rr.size()};
    for (size_t s = 0; s < totalSize; s += m_maxSeqLen)
    {
        size_t end{std::min(s + m_maxSeqLen, totalSize)};
        FloatVector vvBatch(vv.begin() + s, vv.begin() + end);
        FloatVector rrBatch(rr.begin() + s, rr.begin() + end);
        FloatVector cosBatch(cosTheta.begin() + s, cosTheta.begin() + end);
        FloatVector sinBatch(sinTheta.begin() + s, sinTheta.begin() + end);
        FloatVector adcBatch(adc.begin() + s, adc.begin() + end);
        FloatVector widthBatch(width.begin() + s, width.begin() + end);

//        vvBatch.resize(m_maxSeqLen, 0.f);
//        rrBatch.resize(m_maxSeqLen, 0.f);
//        cosBatch.resize(m_maxSeqLen, 0.f);
//        sinBatch.resize(m_maxSeqLen, 0.f);
//        adcBatch.resize(m_maxSeqLen, 0.f);
//        widthBatch.resize(m_maxSeqLen, 0.f);

        long batchSize = static_cast<long>(rrBatch.size());
        std::vector<LArDLHelper::TorchInput> tensors{
            torch::from_blob(vvBatch.data(), {batchSize}),
            torch::from_blob(rrBatch.data(), {batchSize}),
            torch::from_blob(cosBatch.data(), {batchSize}),
            torch::from_blob(sinBatch.data(), {batchSize}),
            torch::from_blob(widthBatch.data(), {batchSize}),
            torch::from_blob(adcBatch.data(), {batchSize})
        };

        // Run the input through the trained model
        LArDLHelper::TorchInputVector inputs;
        LArDLHelper::TorchInput input{torch::stack(tensors, 1)};
        input = input.unsqueeze(0); // Add batch dimension
        std::cout << "Input size: " << input.sizes() << std::endl;
        inputs.push_back(input);
        LArDLHelper::TorchDict output{LArDLHelper::CreateOutputDict()};

        auto mem_func = []()
        {
            std::ifstream statm("/proc/self/statm");
            size_t size, resident;
            statm >> size >> resident;
            return resident * sysconf(_SC_PAGESIZE) / (1024 * 1024);
        };
        size_t before{mem_func()};
        auto t0{std::chrono::high_resolution_clock::now()};
        LArDLHelper::Forward(m_model, inputs, output);
        auto t1{std::chrono::high_resolution_clock::now()};
        size_t after{mem_func()};
        std::cout << "Memory usage: " << (after - before) << " MB with peak: " << after << " MB" << std::endl;
        double duration{std::chrono::duration<double, std::milli>(t1 - t0).count()};
        std::cout << "DL inference time: " << duration << " ms for " << batchSize << " hits" << std::endl;

        // Get the outputs as tensors - we're assuming batch size of 1, hence squeeze(0) for convenience
        LArDLHelper::TorchOutput betaLogits{output.at("beta").toTensor().squeeze(0)}; // [N, 1]
        LArDLHelper::TorchOutput embed{output.at("embed").toTensor().squeeze(0)}; // [N, D]
        LArDLHelper::TorchOutput classLogits{output.at("class_logits").toTensor().squeeze(0)}; // [N, num_classes]
        std::cout << "Output beta size: " << betaLogits.sizes() << std::endl;

        // Convert logits to probabilities
        LArDLHelper::TorchOutput classProbs{torch::softmax(classLogits, 1)};
        // Combine classes 0 and 1 into a single class, in place, representing track-like hits
        classProbs.index({torch::indexing::Slice(), 0}) += classProbs.index({torch::indexing::Slice(), 1});
        classProbs.index({torch::indexing::Slice(), 1}) = 0;
        LArDLHelper::TorchOutput classes{classProbs.argmax(1)}; // [num_cps]
        // If argmax is 3, then check which is large of 0 and 2, and switch to that
        for (long i = 0; i < classes.size(0); ++i)
        {
            if (classes[i].item<int>() == 3)
            {
                const float trackProb{classProbs[i][0].item<float>()};
                const float showerProb{classProbs[i][2].item<float>()};
                if (trackProb > showerProb)
                    classes[i] = 0;
                else
                    classes[i] = 2;
            }
        }

        std::map<int, CaloHitList> classHitsMap{{MIP, CaloHitList()}, {HIP, CaloHitList()}, {SHOWER, CaloHitList()},
            {LOW_E, CaloHitList()}};
        for (long i = 0; i < classes.size(0); ++i)
        {
            const int predictedClass{classes[i].item<int>()};
            const CaloHit *const pCaloHit{caloHitVector[s + i]};
            classHitsMap[predictedClass].push_back(pCaloHit);
            LArCaloHit *const pLArCaloHit{dynamic_cast<LArCaloHit *>(const_cast<CaloHit *>(pCaloHit))};
            if (pLArCaloHit)
            {
                const float denominator{classProbs[i][0].item<float>() + classProbs[i][2].item<float>()};
                const float trackProb{denominator > 0.f ? classProbs[i][0].item<float>() / denominator : 1.f};
                const float showerProb{denominator > 0.f ? classProbs[i][2].item<float>() / denominator : 0.f};
                pLArCaloHit->SetTrackProbability(trackProb);
                pLArCaloHit->SetShowerProbability(showerProb);
            }
        }

        for (const auto &[predictedClass, caloHitList] : classHitsMap)
        {
            std::string categoryName;
            switch (predictedClass + 1)
            {
                case MIP:
                    categoryName = "MIP";
                    break;
                case HIP:
                    categoryName = "HIP";
                    break;
                case SHOWER:
                    categoryName = "SHOWER";
                    break;
                case LOW_E:
                    categoryName = "LOW_E";
                    break;
                default:
                    break;
            }
            // Visualise calo hits for class
            if (!caloHitList.empty())
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, categoryName,
                    predictedClass + 1 == MIP ? BLUE : RED));
            }
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::FoldToLeading(const lar_content::LArMCParticleHelper::MCContributionMap &mcHitMap, const MCFoldingMap &mcFoldingMap,
    lar_content::LArMCParticleHelper::MCContributionMap &leadingHitMap) const
{
    leadingHitMap.clear();
    for (const auto &[pMC, caloHitList] : mcHitMap)
    {
        const MCParticle *const pLeading{mcFoldingMap.at(pMC)};
        leadingHitMap[pLeading].insert(leadingHitMap[pLeading].end(), caloHitList.begin(), caloHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::ConsolidateInstances(const lar_content::LArMCParticleHelper::MCContributionMap &leadingHitMap, 
    lar_content::LArMCParticleHelper::MCContributionMap &instanceHitMap) const
{
    instanceHitMap.clear();
    for (const auto &[pMC, caloHitList] : leadingHitMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC)
        {
            // Ensure delta ray hits along (colinear) muon/pion tracks are assigned to the parent particle
            const int pdg{std::abs(pLArMC->GetParticleId())};
            if (pdg == E_MINUS || pdg == PHOTON)
            {
                CaloHitList particleOwnedHits, parentOwnedHits;
                const LArMCParticle *const pParent{dynamic_cast<const LArMCParticle *>(this->AllocateHitOwner(pLArMC, leadingHitMap,
                    particleOwnedHits, parentOwnedHits))};
                if (!parentOwnedHits.empty())
                    instanceHitMap[pParent].insert(instanceHitMap[pParent].end(), parentOwnedHits.begin(), parentOwnedHits.end());
                if (!particleOwnedHits.empty())
                    instanceHitMap[pLArMC].insert(instanceHitMap[pLArMC].end(), particleOwnedHits.begin(), particleOwnedHits.end());
            }
            else
            {
                if (!caloHitList.empty())
                    instanceHitMap[pLArMC].insert(instanceHitMap[pLArMC].end(), caloHitList.begin(), caloHitList.end());
            }
        }
    }
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

    // Here we will consider Michel electrons to be those that come from muon decays, pion or kaon decays, given the topological similarity
    return (parentPdg == MU_MINUS || parentPdg == PI_PLUS || parentPdg == K_PLUS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlHitTrackShowerIdAlgorithm::IsShower(const LArMCParticle *const pLArMC) const
{
    if (!pLArMC)
        return false;

    const MCProcess process{pLArMC->GetProcess()};
    const int pdg{std::abs(pLArMC->GetParticleId())};

    // Account for the rare(ish) pi zero and K0 long decays involving electrons
    const MCParticleList &parentList{pLArMC->GetParentList()};
    if (!parentList.empty())
    {
        const int parentPdg{std::abs(parentList.front()->GetParticleId())};
        if (((parentPdg == PI_ZERO) || (parentPdg == K_LONG)) && pdg == E_MINUS && process == MC_PROC_DECAY)
            return true;
    }

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

const LArMCParticle *DlHitTrackShowerIdAlgorithm::AllocateHitOwner(const pandora::MCParticle *const pParticle,
    const LArMCParticleHelper::MCContributionMap &mcHitMap, CaloHitList &particleOwnedHits, CaloHitList &parentOwnedHits) const
{
    const MCParticleList &parentList{pParticle->GetParentList()};
    if (parentList.empty())
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    const LArMCParticle *pParent{dynamic_cast<const LArMCParticle *>(parentList.front())};
    if (!(this->IsHip(pParent) || this->IsMip(pParent)))
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
    const CaloHitList &particleHits{mcHitMap.at(pParticle)};
    // Not guaranteed that a parent particle has hits, and so may not be in the hit map
    if (mcHitMap.find(pParent) == mcHitMap.end())
    {
        particleOwnedHits = mcHitMap.at(pParticle);
        return nullptr;
    }
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

StatusCode DlHitTrackShowerIdAlgorithm::GetVertexPlanarCoordinates(const HitType view, float &vx, float &vz) const
{
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_vertexListName, pVertexList));
    if (pVertexList->empty())
    {
        std::cout << "DlHitTrackShowerIdAlgorithm: No vertices found in the vertex list" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }
    const CartesianVector &vertex{m_vertexRelative ? pVertexList->front()->GetPosition() : CartesianVector(0.f, 0.f, 0.f)};
    vx = vertex.GetX();
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    switch (view)
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

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::GetCoordinateExtrema(const CaloHitList &caloHitList, const float vx, const float vz, Bounds &bounds) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        bounds.xMin = std::min(pos.GetX(), bounds.xMin);
        bounds.xMax = std::max(pos.GetX(), bounds.xMax);
        bounds.zMin = std::min(pos.GetZ(), bounds.zMin);
        bounds.zMax = std::max(pos.GetZ(), bounds.zMax);
        if (m_vertexRelative)
        {
            const float x{pos.GetX() - vx};
            const float z{pos.GetZ() - vz};
            const float r{std::sqrt(x * x + z * z)};
            bounds.rMax = std::max(r, bounds.rMax);
        }
    }
    bounds.xRange = bounds.xMax - bounds.xMin > 0 ? bounds.xMax - bounds.xMin : 1.f;
    bounds.zRange = bounds.zMax - bounds.zMin > 0 ? bounds.zMax - bounds.zMin : 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlHitTrackShowerIdAlgorithm::PopulateInputVectors(const CaloHitList &caloHitList, const Bounds &bounds, const float vx, const float vz,
    FloatVector &xx, FloatVector &zz, FloatVector &rr, FloatVector &cosTheta, FloatVector &sinTheta, FloatVector &vv, FloatVector &adc,
    FloatVector &width, CaloHitVector &sortedCaloHitList) const
{
    if (caloHitList.empty())
        return;

    const size_t nHits{caloHitList.size()};

    if (!m_polarCoords)
    {
        xx.resize(nHits, 0);
        zz.resize(nHits, 0);
    }
    else
    {
        rr.resize(nHits, 0);
        cosTheta.resize(nHits, 0);
        sinTheta.resize(nHits, 0);
    }
    vv.resize(nHits, 0);
    adc.resize(nHits, 0);
    width.resize(nHits, 0);

    size_t i{0};
    for (const CaloHit *const pCaloHit : caloHitList)
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
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }
        if (!m_polarCoords)
        {
            xx[i] = (x - bounds.xMin) / bounds.xRange;
            zz[i] = (z - bounds.zMin) / bounds.zRange;
        }
        else
        {
            const float dx{x - vx};
            const float dz{z - vz};
            const float r{std::sqrt(dx * dx + dz * dz)};
            if (r > 0.f)
            {
                rr[i] = r / bounds.rMax;
                // Scale to [0,1]
                cosTheta[i] = 0.5f * (1.f + dx / r);
                sinTheta[i] = 0.5f * (1.f + dz / r);
            }
            else
            {
                rr[i] = 0.f;
                // Scale to [0,1]
                cosTheta[i] = 0.5f * (1.f + 1.f);
                sinTheta[i] = 0.5f * (1.f + 0.f);
            }
        }
        adc[i] = std::log(1.f + std::min(pCaloHit->GetInputEnergy(), m_adcPeak * m_maxAdcFactor)) / std::log(1.f + m_adcPeak * m_maxAdcFactor);
        width[i] = pCaloHit->GetCellSize1() / bounds.xRange;
        ++i;
    }

    if (m_polarCoords)
    {
        // Sort the vectors by r
        IndexVector indices(adc.size());
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

        sortedCaloHitList.clear();
        sortedCaloHitList.insert(sortedCaloHitList.end(), caloHitList.begin(), caloHitList.end());
        reorder(sortedCaloHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlHitTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
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
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_model);
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AdcPeak", m_adcPeak));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAdcFactor", m_maxAdcFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
