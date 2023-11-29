/**
 *  @file   larpandoradlcontent/LArThreeDReco/DlPfoEvaluatorAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArThreeDReco/DlPfoEvaluatorAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlPfoEvaluatorAlgorithm::DlPfoEvaluatorAlgorithm() :
    m_caloHitListName{"CaloHitList2D"},
    m_purityThreshold{0.f},
    m_completenessThreshold{0.f},
    m_useAdcWeighting{false},
    m_trainingMode{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlPfoEvaluatorAlgorithm::~DlPfoEvaluatorAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreePrefix + "_u", m_rootFileName, "RECREATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreePrefix + "_v", m_rootFileName, "UPDATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreePrefix + "_w", m_rootFileName, "UPDATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlPfoEvaluatorAlgorithm::~DlPfoEvaluatorAlgorithm - Unable to write to ROOT tree" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoEvaluatorAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoEvaluatorAlgorithm::PrepareTrainingSample()
{
    MapOfMCHitMaps mcToAllHitsMap;
    const bool usePurity{m_purityThreshold > 0.f};
    const float threshold{std::max(m_purityThreshold, m_completenessThreshold)};
    if (usePurity)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
        this->GetMCHitsMap(*pCaloHitList, mcToAllHitsMap);
    }

    for (const std::string &listname : m_pfoListNames)
    {
        const PfoList *pPfoList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pPfoList));

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            if (usePurity)
            {
                this->CreateTrainingExample(pPfo, TPC_VIEW_U, this->GetPurity(pPfo, TPC_VIEW_U), threshold);
                this->CreateTrainingExample(pPfo, TPC_VIEW_V, this->GetPurity(pPfo, TPC_VIEW_V), threshold);
                this->CreateTrainingExample(pPfo, TPC_VIEW_W, this->GetPurity(pPfo, TPC_VIEW_W), threshold);
            }
            else
            {
                this->CreateTrainingExample(pPfo, TPC_VIEW_U, this->GetCompleteness(pPfo, TPC_VIEW_U, mcToAllHitsMap), threshold);
                this->CreateTrainingExample(pPfo, TPC_VIEW_V, this->GetCompleteness(pPfo, TPC_VIEW_V, mcToAllHitsMap), threshold);
                this->CreateTrainingExample(pPfo, TPC_VIEW_W, this->GetCompleteness(pPfo, TPC_VIEW_W, mcToAllHitsMap), threshold);
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoEvaluatorAlgorithm::Infer()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DlPfoEvaluatorAlgorithm::GetPurity(const ParticleFlowObject *const pPfo, const HitType view)
{
    CaloHitList pfoHitList;
    LArPfoHelper::GetCaloHits(pPfo, view, pfoHitList);
    float weight{0.f};
    const MCParticle *const pMainMC{this->GetMainMCParticle(pfoHitList, weight)};

    if (pMainMC && weight > 0.f)
        return weight * this->GetFactor(pfoHitList);
    else
        return -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DlPfoEvaluatorAlgorithm::GetCompleteness(const ParticleFlowObject *const pPfo, const HitType view, const MapOfMCHitMaps &mcToAllHitsMap)
{
    CaloHitList pfoHitList;
    LArPfoHelper::GetCaloHits(pPfo, view, pfoHitList);
    float weight{0.f};
    const MCParticle *const pMainMC{this->GetMainMCParticle(pfoHitList, weight)};

    if (pMainMC && weight > 0.f)
        return weight * this->GetFactor(mcToAllHitsMap.at(view).at(pMainMC));
    else
        return -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DlPfoEvaluatorAlgorithm::GetFactor(const CaloHitList &caloHitList)
{
    float denomenator{0.f};
    for (const CaloHit *const pCaloHit : caloHitList)
        denomenator += m_useAdcWeighting ? pCaloHit->GetInputEnergy() : 1.f;
    // ATTN - weight > 0 check in calling functions means we don't need to check for zero denomenator here
    return 1.f / denomenator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *DlPfoEvaluatorAlgorithm::GetMainMCParticle(const CaloHitList &pfoHitList, float &weight)
{
    MCHitMap mcToPfoHitsMap;
    this->GetMCHitsMap(pfoHitList, mcToPfoHitsMap);
    const MCParticle *pMainMC{nullptr};
    for (auto &[pMC, sharedHitList] : mcToPfoHitsMap)
    {
        float totalWeight{0.f};
        for (const CaloHit *const pCaloHit : sharedHitList)
            totalWeight += m_useAdcWeighting ? pCaloHit->GetInputEnergy() : 1.f;
        if (totalWeight > weight)
        {
            weight = totalWeight;
            pMainMC = pMC;
        }
    }

    return pMainMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoEvaluatorAlgorithm::CreateTrainingExample(const ParticleFlowObject *const pPfo, const HitType view, const float metric, const float threshold)
{
    if (metric < 0)
        return;
    FloatVector x, z, q;
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, view, caloHitList);
    CaloHitVector caloHitVector;
    for (const CaloHit *const pCaloHit : caloHitList)
        caloHitVector.emplace_back(pCaloHit);
    const Vertex *const pVertex{LArPfoHelper::GetVertex(pPfo)};
    if (!pVertex)
        return;
    const CartesianVector &vertex3D{pVertex->GetPosition()};
    const CartesianVector &vertex2D{LArGeometryHelper::ProjectPosition(this->GetPandora(), vertex3D, view)};

    auto sortByDistance = [&](const CaloHit *const pCaloHit1, const CaloHit *const pCaloHit2)
    {
        const float d1{pCaloHit1->GetPositionVector().GetDistanceSquared(vertex2D)};
        const float d2{pCaloHit2->GetPositionVector().GetDistanceSquared(vertex2D)};
        return d1 < d2;
    };

    std::sort(caloHitVector.begin(), caloHitVector.end(), sortByDistance);
    for (const CaloHit *const pCaloHit : caloHitVector)
    {
        // Should sort these according to some distance metric
        x.emplace_back(pCaloHit->GetPositionVector().GetX());
        z.emplace_back(pCaloHit->GetPositionVector().GetZ());
        q.emplace_back(pCaloHit->GetMipEquivalentEnergy());
    }

    const std::string treeName{m_rootTreePrefix + (view == TPC_VIEW_U ? "_u" : view == TPC_VIEW_V ? "_v" : "_w")};
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "x", x));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "z", z));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "q", q));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), treeName, "is_good", metric >= threshold));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), treeName));
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoEvaluatorAlgorithm::GetMCHitsMap(const CaloHitList caloHitList, MCHitMap &mcHitsMap)
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            mcHitsMap[pMC].emplace_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoEvaluatorAlgorithm::GetMCHitsMap(const CaloHitList caloHitList, MapOfMCHitMaps &mcHitsMap)
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            mcHitsMap[pCaloHit->GetHitType()][pMC].emplace_back(pCaloHit);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoEvaluatorAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreePrefix", m_rootTreePrefix));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName",
            m_caloHitListName));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PurityThreshold",
            m_purityThreshold));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CompletenessThreshold",
            m_completenessThreshold));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseAdcWeighting",
            m_useAdcWeighting));

        if (m_purityThreshold == 0.f && m_completenessThreshold == 0.f)
        {
            std::cout << "DlPfoEvaluatorAlgorithm::ReadSettings - One of purity or completeness must have a non-zero threshold" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
        else if (m_purityThreshold > 0.f && m_completenessThreshold > 0.f)
        {
            std::cout << "DlPfoEvaluatorAlgorithm::ReadSettings - Only one of purity or completeness should have a non-zero threshold" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
