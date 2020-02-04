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

#include <chrono>

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

bool isNeutronDaughter(const MCParticle *const mcParticle)
{
    const MCParticle* particle = mcParticle;
    while(!particle->GetParentList().empty())
    {   
        if(particle->GetParentList().size() > 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const MCParticle* pParent = *(particle->GetParentList().begin());
        const int pdg = std::abs(pParent->GetParticleId());
        if (pdg == 2112)
        {
            return true;
        }   
        particle = pParent;
    }

    return false;
}


StatusCode DeepLearningTrackShowerIdAlgorithm::Run()
{
    if (m_useTrainingMode)
        return Train();
    else
        return Infer();
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

StatusCode DeepLearningTrackShowerIdAlgorithm::Train()
{
    for (const std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW) return STATUS_CODE_NOT_ALLOWED;

        std::string trainingOutputFileName(m_trainingOutputFile);
        LArMvaHelper::MvaFeatureVector featureVector;

        if (isU) trainingOutputFileName += "_CaloHitListU.txt";
        else if (isV) trainingOutputFileName += "_CaloHitListV.txt";
        else if (isW) trainingOutputFileName += "_CaloHitListW.txt";

        featureVector.push_back(static_cast<double>(pCaloHitList->size()));

        int discardedHits = 0;
        LArMCParticleHelper::PrimaryParameters parameters;
        // Only care about reconstructability with respect to the current view, so skip good view check
        parameters.m_minHitsForGoodView = 0;
        // Turn off max photo propagation for now, only care about killing off daughters of neutrons
        parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList,
                parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            int isReconstructable = 1;
            int pdg(-1);

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                // Throw away non-reconstructable hits
                if (targetMCParticleToHitsMap.find(pMCParticle) == targetMCParticleToHitsMap.end())
                {
                    ++discardedHits;
                    isReconstructable = 0;
                }
                pdg = pMCParticle->GetParticleId();
                if (isReconstructable)
                {
                    if(isNeutronDaughter(pMCParticle))
                    {
                        isReconstructable = 0;
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
            featureVector.push_back(static_cast<double>(isReconstructable));
        }

        std::cout << "Discarding " << discardedHits << " of " << pCaloHitList->size() <<  " hits as non-reconstructable" << std::endl;

        LArMvaHelper::ProduceTrainingExample(trainingOutputFileName, true, featureVector);
    }
    return STATUS_CODE_SUCCESS;
}

StatusCode DeepLearningTrackShowerIdAlgorithm::Infer()
{
    std::cout << "system_clock " << std::chrono::system_clock::period::num << " " <<
        std::chrono::system_clock::period::den << " steady: " << std::boolalpha <<
        std::chrono::system_clock::is_steady << std::endl;

    std::cout << "high_resolution_clock " << std::chrono::high_resolution_clock::period::num << " " <<
        std::chrono::high_resolution_clock::period::den << " steady: " << std::boolalpha <<
        std::chrono::high_resolution_clock::is_steady << std::endl;

    std::cout << "steady_clock " << std::chrono::steady_clock::period::num << " " <<
        std::chrono::steady_clock::period::den << " steady: " << std::boolalpha <<
        std::chrono::steady_clock::is_steady << std::endl;
  
    auto start = std::chrono::steady_clock::now();
    std::shared_ptr<torch::jit::script::Module> pModule(nullptr);
    try
    {
        pModule = torch::jit::load(m_modelFileName);
    }
    catch (...)
    {
        std::cout << "Error loading the PyTorch module" << std::endl;
        return STATUS_CODE_FAILURE;
    }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << "Network setup time: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;

    if (m_visualize)
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    for (const std::string listName : m_caloHitListNames)
    {
        start = std::chrono::steady_clock::now();
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        const bool isU(pCaloHitList->front()->GetHitType() == TPC_VIEW_U ? true : false);
        const bool isV(pCaloHitList->front()->GetHitType() == TPC_VIEW_V ? true : false);
        const bool isW(pCaloHitList->front()->GetHitType() == TPC_VIEW_W ? true : false);

        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

        const float tileSize = 128.f;
        const int imageWidth = 256;
        const int imageHeight = 256;
        // Get bounds of hit region
        float xMin = std::numeric_limits<float>::max(); float xMax = 0.f;
        float zMin = std::numeric_limits<float>::max(); float zMax = 0.f;
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());
            if (x < xMin) xMin = x;
            if (x > xMax) xMax = x;
            if (z < zMin) zMin = z;
            if (z > zMax) zMax = z;
        }
        const float xRange = xMax - xMin;
        const float zRange = zMax - zMin;
        int nTilesX = static_cast<int>(std::ceil(xRange / tileSize));
        int nTilesZ = static_cast<int>(std::ceil(zRange / tileSize));
        // Need to add 1 to number of tiles for ranges exactly matching tile size and for zero ranges
        if (std::fmod(xRange, tileSize) == 0.f) ++nTilesX;
        if (std::fmod(zRange, tileSize) == 0.f) ++nTilesZ;

        const int nTiles = nTilesX * nTilesZ;

        // Normalisation value from training set - should store and load this from file
        const float PIXEL_VALUE = (255.0 - 0.558497428894043) / 11.920382499694824;
        // Populate tensor - may need to split the batch into separate single image batches
        torch::Tensor input = torch::zeros({nTiles, 1, imageHeight, imageWidth});
        typedef std::map<const CaloHit*, std::tuple<int, int, int, int>> CaloHitToPixelMap;
        CaloHitToPixelMap caloHitToPixelMap;
        auto accessor = input.accessor<float, 4>();
        end = std::chrono::steady_clock::now();
        diff = end - start;
        std::cout << "#Hits: " << pCaloHitList->size() << " #TilesX: " <<
            nTilesX << " #TilesZ: " << nTilesZ << " #Tiles: " << nTiles << std::endl;
        std::cout << "Prep time: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;
        start = std::chrono::steady_clock::now();
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());
            // Determine which tile the hit will be assigned to
            const int tileX = static_cast<int>(std::floor((x - xMin) / tileSize));
            const int tileZ = static_cast<int>(std::floor((z - zMin) / tileSize));
            const int tile = tileZ * nTilesX + tileX;
            // Determine hit position within the tile
            const float localX = std::fmod(x - xMin, tileSize);
            const float localZ = std::fmod(z - zMin, tileSize);
            // Determine hit pixel within the tile
            const int pixelX = static_cast<int>(std::floor(localX * imageHeight / tileSize));
            const int pixelZ = static_cast<int>(std::floor(localZ * imageWidth / tileSize));
            accessor[tile][0][pixelX][pixelZ] = PIXEL_VALUE;
            caloHitToPixelMap.insert(std::make_pair(pCaloHit, std::make_tuple(tileX, tileZ, pixelX, pixelZ)));
        }

        // Pass as input the input Tensor containing the calo hit picture
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(input);
        end = std::chrono::steady_clock::now();
        diff = end - start;
        std::cout << "Image processing time: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;

        // Run the input through the trained model and get the output accessor
        // There should be one iteration of this for each image associated with a given tile
        // Need to ensure that each calohit is mapped to the correct tile and extract the
        // track/shower probability from the appropriate pixel within that particular ouput
        // This may need multiple runs of forward on the different tile inputs - essentially
        // this would need to know about the value in the first dimension of the outputAccessor -
        // i.e. the number within the batch
        start = std::chrono::steady_clock::now();
        at::Tensor output = pModule->forward(inputs).toTensor();
        auto outputAccessor = output.accessor<float, 4>();
        end = std::chrono::steady_clock::now();
        diff = end - start;
        std::cout << "Network time: " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;

        start = std::chrono::steady_clock::now();
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            auto pixelMap = caloHitToPixelMap.at(pCaloHit);
            const int tileX(std::get<0>(pixelMap));
            const int tileZ(std::get<1>(pixelMap));
            const int tile = tileZ * nTilesX + tileX;
            const int pixelX(std::get<2>(pixelMap));
            const int pixelZ(std::get<3>(pixelMap));

            object_creation::CaloHit::Metadata metadata;
            metadata.m_propertiesToAdd["Pshower"] = outputAccessor[tile][1][pixelX][pixelZ];
            metadata.m_propertiesToAdd["Ptrack"] = outputAccessor[tile][2][pixelX][pixelZ];
            const StatusCode &statusCode(PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
            if (statusCode != STATUS_CODE_SUCCESS)
                std::cout << "Cannot set calo hit meta data" << std::endl;
        }
        end = std::chrono::steady_clock::now();
        diff = end - start;
        std::cout << "Readout time " << std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;

        /*

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
        }*/
    }

    if (m_visualize)
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
