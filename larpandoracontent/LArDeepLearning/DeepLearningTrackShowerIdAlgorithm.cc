/**
 *  @file   larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower id algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include <torch/script.h>

#include "larpandoracontent/LArDeepLearning/DeepLearningTrackShowerIdAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

DeepLearningTrackShowerIdAlgorithm::DeepLearningTrackShowerIdAlgorithm() :
    m_xMin(-420),
    m_zMinU(-350),
    m_zMinV(0),
    m_zMinW(-25),
    m_nBins(512),
    m_visualize(false),
    m_useTrainingMode(false),
    m_trainingOutputFile("")
{
    const float span(980);
    m_xMax = m_xMin + span;
    m_zMaxU = m_zMinU + span;
    m_zMaxV = m_zMinV + span;
    m_zMaxW = m_zMinW + span;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool isIsolatedShortTrack(const MCParticle *const mcParticle)
{
    const MCParticle* particle = mcParticle;
    while(!particle->GetParentList().empty())
    {   
        if(particle->GetParentList().size() > 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const MCParticle* pParent = *(particle->GetParentList().begin());
        const int pdg = std::abs(pParent->GetParticleId());
        if(pdg == 211 || pdg == 2112)
        {
            CartesianVector displacement = mcParticle->GetEndpoint() - mcParticle->GetVertex();
            if(displacement.GetMagnitude() < 10)
                return true;
        }   
        particle = pParent;
    }

    return false;
}


StatusCode DeepLearningTrackShowerIdAlgorithm::Run()
{
    if (m_useTrainingMode)
    {
        for (const std::string listName : m_caloHitListNames)
        {
            const CaloHitList *pCaloHitList(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

            const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
            const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
            const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

            if (!isU && !isV && !isW)
                return STATUS_CODE_NOT_ALLOWED;

            std::string trainingOutputFileName(m_trainingOutputFile);
            LArMvaHelper::MvaFeatureVector featureVector;

            if (isU) trainingOutputFileName += "_CaloHitListU.txt";
            else if (isV) trainingOutputFileName += "_CaloHitListV.txt";
            else if (isW) trainingOutputFileName += "_CaloHitListW.txt";

            featureVector.push_back(static_cast<double>(pCaloHitList->size()));

            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                int nuanceCode(3000);
                int pdg(-1);

                try
                {
                    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                    nuanceCode = LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
                    pdg = pMCParticle->GetParticleId();
                    if(std::abs(pdg) != 11 && std::abs(pdg) != 22)
                    {   // Lots of neutrons and some pions bumping into protons producing isolated short tracks - force them to be showers
                        if(isIsolatedShortTrack(pMCParticle))
                        {
                           pdg = 11;
                        }
                    }
                }
                catch (...)
                {
                    continue;
                }

                featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
                featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetY()));
                featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
                featureVector.push_back(static_cast<double>(pdg));
                featureVector.push_back(static_cast<double>(nuanceCode));
            }

            LArMvaHelper::ProduceTrainingExample(trainingOutputFileName, true, featureVector);
        }
        return STATUS_CODE_SUCCESS;
    }

    // Load the model.pt file.
    std::shared_ptr<torch::jit::script::Module> pModule(nullptr);

    try
    {
        pModule = torch::jit::load(m_modelFileName);
    }
    catch (const c10::Error &e)
    {
        std::cout << "Error loading the PyTorch module" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    if (m_visualize)
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        const float zMin(isU ? m_zMinU : (isV ? m_zMinV : m_zMinW));
        const float zMax(isU ? m_zMaxU : (isV ? m_zMaxV : m_zMaxW));

        const float xSpan(m_xMax - m_xMin), zSpan(zMax - zMin);

        typedef std::map<const CaloHit*, std::pair<int, int>> CaloHitToBinMap;
        CaloHitToBinMap caloHitToBinMap;

        // Start with RGB picture of black pixels.  Four indices: first default size 1, second index is RGB indices, third is xBin,
        // fourth is zBin
        torch::Tensor input = torch::zeros({1, 3, m_nBins, m_nBins});
        auto accessor = input.accessor<float, 4>();

        // Create a map of calo hits to x/z bin values.  Set the output track shower id of the pixel using the RGB values at the pixel
        // containing the calo hit
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());

            const int xBin(std::floor((x-m_xMin)*m_nBins/xSpan));
            const int zBin(std::floor((z-zMin)*m_nBins/zSpan));

            // ATTN: Set pixels containing a calo hit to white
            if (xBin >= 0 && xBin <= m_nBins && zBin >= 0 && zBin <= m_nBins)
            {
                caloHitToBinMap.insert(std::make_pair(pCaloHit, std::make_pair(xBin, zBin)));
                accessor[0][0][xBin][zBin] = 1;
                accessor[0][1][xBin][zBin] = 1;
                accessor[0][2][xBin][zBin] = 1;
            }
        }

        // Pass as input the input Tensor containing the calo hit picture
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(input);

        // Run the input through the trained model and get the output accessor
        at::Tensor output = pModule->forward(inputs).toTensor();
        auto outputAccessor = output.accessor<float, 4>();

        // Colour in the shower and track bits (and other) in a visual display for first performance inspection
        CaloHitList showerHits, trackHits, other;

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            if (caloHitToBinMap.find(pCaloHit) == caloHitToBinMap.end())
            {
                other.push_back(pCaloHit);
                continue;
            }

            const int xBin(caloHitToBinMap.at(pCaloHit).first);
            const int zBin(caloHitToBinMap.at(pCaloHit).second);

            // Is the R value bigger than the B value.  In training the target picture was coloured such that showers were red and tracks blue
            const bool isShower(outputAccessor[0][0][xBin][zBin] > outputAccessor[0][2][xBin][zBin] ? true : false);
            object_creation::CaloHit::Metadata metadata;

            if (isShower)
            {
                metadata.m_propertiesToAdd["IsShower"] = 1.f;
                showerHits.push_back(pCaloHit);
            }
            else
            {
                metadata.m_propertiesToAdd["IsTrack"] = 1.f;
                trackHits.push_back(pCaloHit);
            }

            const StatusCode &statusCode(PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));

            if (statusCode != STATUS_CODE_SUCCESS)
                std::cout << "Cannot set calo hit meta data" << std::endl;
        }

        if (m_visualize)
        {
            const std::string trackListName("TrackHits_" + listName);
            const std::string showerListName("ShowerHits_" + listName);
            const std::string otherListName("OtherHits_" + listName);
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHits, trackListName, BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHits, showerListName, RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &other, otherListName, BLACK));
        }
    }

    if (m_visualize)
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", m_modelFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Visualize", m_visualize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageXMin", m_xMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageXMax", m_xMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMinU", m_zMinU));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMaxU", m_zMaxU));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMinV", m_zMinV));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMaxV", m_zMaxV));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMinW", m_zMinW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ImageZMaxW", m_zMaxW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumberOfBins", m_nBins));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
