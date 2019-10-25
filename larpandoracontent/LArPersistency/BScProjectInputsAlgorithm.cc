/**
 *  @file   PandoraSDK/src/Persistency/BScProjectInputsAlgorithm.cc
 *
 *  @brief  Implementation of the event writing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArPersistency/BScProjectInputsAlgorithm.h"

#include "Pandora/PandoraInternal.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "Objects/OrderedCaloHitList.h"

#include <fstream>

using namespace pandora;
using namespace lar_content;

namespace lar_content
{

BScProjectInputsAlgorithm::BScProjectInputsAlgorithm():
    m_inputCaloHitListName(""),
    m_outputFilename("")
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

BScProjectInputsAlgorithm::~BScProjectInputsAlgorithm()
{
    if(m_file.is_open())
    {
        m_file.close();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BScProjectInputsAlgorithm::Run()
{
/*    const CaloHitList *pCaloHitList(nullptr);
    std::cout << "Name: " << m_inputCaloHitListName << std::endl;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(
        *this, m_inputCaloHitListName, pCaloHitList));*/


    ///////////////
    
    if(!m_file.is_open())
    {
        return STATUS_CODE_FAILURE;
    }

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(
        *this, "NeutrinoParticles3D", pPfoList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(
        *this, m_inputCaloHitListName, pCaloHitList));

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    // Mapping target MCParticles -> truth associated Hits
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList,
        pCaloHitList, LArMCParticleHelper::PrimaryParameters(),
        LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList,
        targetMCParticleToHitsMap, pfoToHitsMap);

    // Last step
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap,
        {targetMCParticleToHitsMap}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);

    if (pPfoList->size() > 1)
    {
        std::cout << "WARNING: Found more than 1 Pfo!" << std::endl;
    }
    for (const Pfo *const pPfo : *pPfoList)
    {
        auto vertices = pPfo->GetVertexList();
        if (vertices.size() == 1)
        {            
            const CartesianVector trueVertex = vertices.front()->GetPosition();
            m_file << trueVertex.GetX() << "," << trueVertex.GetY() << "," <<
                trueVertex.GetZ() << ",";
        }
        else
        {
            std::cout << "WARNING: Multiple vertices found" << std::endl;
            m_file << "," << "," << ",";
        }
        
        MCParticleList mcPrimaryList;
        for (const auto &mapEntry : mcToPfoHitSharingMap)
        {
            mcPrimaryList.push_back(mapEntry.first);
        }
        mcPrimaryList.sort(LArMCParticleHelper::SortByMomentum);
        const LArInteractionTypeHelper::InteractionType interactionType(
            LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
        std::string int_type = LArInteractionTypeHelper::ToString(interactionType);
        std::cout << int_type << std::endl;
        m_file << int_type << ",";

        CaloHitList wHitsInPfo;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, wHitsInPfo);
        
        const LArMCParticleHelper::MCParticleToSharedHitsVector
            &mcParticleToSharedHitsVector(pfoToMCHitSharingMap.at(pPfo));
        CaloHitList associatedMCHitsW;
        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
        {
            const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);
            
            for (const CaloHit *const pCaloHit : associatedMCHits)
            {
                if (TPC_VIEW_W == pCaloHit->GetHitType())
                    associatedMCHitsW.push_back(pCaloHit);
            }            
        }
        
        // Record the total number of calo hits
        m_file << associatedMCHitsW.size() << ",";
        // Loop over all of the W calo hits and store their locations
        for (const CaloHit *const pCaloHit : associatedMCHitsW)
        {
            const CartesianVector hitPos = pCaloHit->GetPositionVector();
            float energy = pCaloHit->GetInputEnergy();
            
            m_file << hitPos.GetX() << " " << hitPos.GetZ() << " " << energy << " ";
        }
        m_file << std::endl;
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BScProjectInputsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND,
        !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName",
        m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND,
        !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName",
        m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND,
        !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_outputFilename));
        
    if(!m_outputFilename.empty())
    {
        m_file.open(m_outputFilename, std::ios::out);
        if(m_file.is_open())
        {
            m_file << "true_vertex_x,true_vertex_y,true_vertex_z,interaction_type,num_hits,hits(x wire energy)" << std::endl;
        }
    }
            
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
