/**
 *  @file   larpandoracontent/LArPersistency/HitWritingAlgorithm.cc
 *
 *  @brief  Implementation of the hit writing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArPersistency/HitWritingAlgorithm.h"

#include "Pandora/PandoraInternal.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "Objects/OrderedCaloHitList.h"

using namespace pandora;

namespace lar_content
{

HitWritingAlgorithm::HitWritingAlgorithm():
    m_inputCaloHitList2DName{"CaloHitList2D"},
    m_outputFilename{""},
    m_treeName{"tree"}
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

HitWritingAlgorithm::~HitWritingAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_outputFilename.c_str(), "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "HitWritingAlgorithm: Unable to write tree to file " << m_outputFilename << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWritingAlgorithm::Run()
{
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(
        *this, m_inputCaloHitList2DName, pCaloHitList2D));
    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(
        *this, pMCParticleList));

    if (!(pCaloHitList2D && pMCParticleList))
        return STATUS_CODE_FAILURE;

    // Select reconstructable hits, folding back to primaries for interaction type determination
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_foldBackHierarchy = true;
    LArMCParticleHelper::MCContributionMap primaryMCParticlesToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, primaryMCParticlesToHitsMap);

    // Determine the interaction type from MC truth
    MCParticleList mcPrimaryList;
    for (const auto [ pMCParticle, hits ] : primaryMCParticlesToHitsMap)
    {   (void) hits;
        mcPrimaryList.push_back(pMCParticle);
    }
    mcPrimaryList.sort(LArMCParticleHelper::SortByMomentum);
    const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    const std::string interactionStr{LArInteractionTypeHelper::ToString(interactionType)};
    const int trueInteraction{static_cast<int>(interactionType)};

    const HitType validViews[]{TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

    if (std::find(std::begin(m_interactions), std::end(m_interactions), interactionStr) != std::end(m_interactions))
    {
        float trueVrtX{0.f}, trueVrtY{0.f}, trueVrtZ{0.f};
        FloatVector xCoordsU, zCoordsU, energiesU;
        FloatVector xCoordsV, zCoordsV, energiesV;
        FloatVector xCoordsW, zCoordsW, energiesW;
        for (const MCParticle *pMCParticle : *pMCParticleList)
        {
            if(LArMCParticleHelper::IsNeutrino(pMCParticle))
            {
                const CartesianVector& trueVertex{pMCParticle->GetVertex()};
                trueVrtX = trueVertex.GetX();
                trueVrtY = trueVertex.GetY();
                trueVrtZ = trueVertex.GetZ();
            }
        }

        for (const std::string &caloHitListName : m_inputCaloHitListNames)
        {
            const CaloHitList *pCaloHitList{nullptr};
            PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, caloHitListName, pCaloHitList));

            if (!pCaloHitList || pCaloHitList->empty())
                continue;

            const HitType hitType((*(pCaloHitList->begin()))->GetHitType());
            if (std::find(std::begin(validViews), std::end(validViews), hitType) == std::end(validViews))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                const CartesianVector hitPos = pCaloHit->GetPositionVector();
                float energy = pCaloHit->GetInputEnergy();

                if (hitType == TPC_VIEW_U)
                {
                    xCoordsU.push_back(hitPos.GetX());
                    zCoordsU.push_back(hitPos.GetZ());
                    energiesU.push_back(energy);
                }
                else if (hitType == TPC_VIEW_V)
                {
                    xCoordsV.push_back(hitPos.GetX());
                    zCoordsV.push_back(hitPos.GetZ());
                    energiesV.push_back(energy);
                }
                else
                {
                    xCoordsW.push_back(hitPos.GetX());
                    zCoordsW.push_back(hitPos.GetZ());
                    energiesW.push_back(energy);
                }
            }
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueInteraction", trueInteraction));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueVrtX", trueVrtX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueVrtY", trueVrtY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueVrtZ", trueVrtZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xCoordsU", &xCoordsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zCoordsU", &zCoordsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energiesU", &energiesU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xCoordsV", &xCoordsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zCoordsV", &zCoordsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energiesV", &energiesV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xCoordsW", &xCoordsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zCoordsW", &zCoordsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energiesW", &energiesW));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWritingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputCaloHitListNames", m_inputCaloHitListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitList2DName", m_inputCaloHitList2DName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_outputFilename));
        
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

