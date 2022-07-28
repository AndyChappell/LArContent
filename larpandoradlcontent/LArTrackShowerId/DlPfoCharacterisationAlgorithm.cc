/**
 *  @file   larpandoracontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoradlcontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.h"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlPfoCharacterisationAlgorithm::DlPfoCharacterisationAlgorithm() :
    m_training{false},
    m_imageWidth{224},
    m_imageHeight{224},
    m_minHitsForGoodView{5}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::Run()
{
    if (m_training)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::Infer()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    PfoList tracksToShowers, showersToTracks;

    const PfoList *pPfoList(nullptr);
    std::map<const std::string, LArDLHelper::TorchModel> models{{"U", m_modelU}, {"V", m_modelV}, {"W", m_modelW}};
    std::map<const std::string, const HitType> views{{"U", HitType::TPC_VIEW_U}, {"V", HitType::TPC_VIEW_V}, {"W", HitType::TPC_VIEW_W}};
    for (const bool inputIsTrack : {true, false})
    {
        const std::string listName{inputIsTrack ? m_trackPfoListName : m_showerPfoListName};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pPfoList));

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            int nTrackVotes{0}, nShowerVotes{0};
            for (const std::string &view : {"U", "V", "W"})
            {
                CaloHitList caloHitList;
                LArPfoHelper::GetCaloHits(pPfo, views[view], caloHitList);
                if (caloHitList.empty())
                    continue;
                LArDLHelper::TorchInput input;
                LArDLHelper::InitialiseInput({1, 3, m_imageHeight, m_imageWidth}, input);
                this->PrepareNetworkInput(caloHitList, input);

                // Run the input through the trained model and get the output accessor
                LArDLHelper::TorchInputVector inputs;
                inputs.emplace_back(input);
                LArDLHelper::TorchOutput output;
                LArDLHelper::Forward(models[view], inputs, output);
                auto outputAccessor = output.accessor<float, 2>();
                // Get class probabilities
                const double outputShower{exp(outputAccessor[0][0])};
                const double outputTrack{exp(outputAccessor[0][1])};
                const double sum{(outputShower + outputTrack) > std::numeric_limits<float>::epsilon() ? outputShower + outputTrack : 1.f};
                const double probShower{outputShower / sum}, probTrack{outputTrack / sum};
                if (probShower > probTrack)
                    ++nShowerVotes;
                else
                    ++nTrackVotes;
                {
                    const std::string netClass{probTrack >= probShower ? "Track" : "Shower"};

                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, view, BLUE));
                    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                }
            }
            std::cout << "Here" << std::endl;
            if (inputIsTrack && nShowerVotes > nTrackVotes)
            {
                this->ChangeCharacterisation(pPfo);
                tracksToShowers.emplace_back(pPfo);
                std::cout << "Changed to shower truth, " << this->GetTrackFraction(pPfo) << std::endl;
            }
            else if (!inputIsTrack && nTrackVotes > nShowerVotes)
            {
                this->ChangeCharacterisation(pPfo);
                showersToTracks.emplace_back(pPfo);
                std::cout << "Changed to track, truth " << this->GetTrackFraction(pPfo) << std::endl;
            }
            else
            {
                std::cout << "Unchanged. Retained as " << inputIsTrack << " truth " << this->GetTrackFraction(pPfo) << std::endl;
            }
        }
    }

    if (!tracksToShowers.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

    if (!showersToTracks.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::PrepareTrainingSample()
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ProcessPfoList(m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ProcessPfoList(m_showerPfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoCharacterisationAlgorithm::PrepareNetworkInput(const CaloHitList &caloHitList, LArDLHelper::TorchInput &input) const
{
    if (caloHitList.empty())
        return;

    const double eps{1.1920929e-7}; // Python float epsilon, used in image padding
    double xMin{0.}, xMax{0.}, zMin{0.}, zMax{0.};

    this->GetHitRegion(caloHitList, xMin, xMax, zMin, zMax);
    double xRange{xMax - xMin}, zRange{zMax - zMin};

    // Enforce a minimum region size
    if (xRange < 112.)
    {
        const double pad{0.5 * (112. - xRange)};
        xMin -= pad;
        xMax += pad;
    }
    if (zRange < 112.)
    {
        const double pad{0.5 * (112. - zRange)};
        zMin -= pad;
        zMax += pad;
    }
    xMin -= eps; xMax += eps;
    zMin -= eps; zMax += eps;
    xRange = xMax - xMin;
    zRange = zMax - zMin;

    // Channel mapping is UVW => RGB (012)
    auto accessor = input.accessor<float, 4>();
    for (const CaloHit *pCaloHit : caloHitList)
    {
        // Determine hit pixel
        const double x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
        const double localX{(x - xMin) / xRange};
        const double localZ{(z - zMin) / zRange};
        const int pixelX{static_cast<int>(std::floor(localX * m_imageWidth))};
        const int pixelZ{static_cast<int>(std::floor(localZ * m_imageHeight))};
        // Repeat channel for RGB input image
        accessor[0][0][pixelX][pixelZ] += pCaloHit->GetInputEnergy();
        accessor[0][1][pixelX][pixelZ] += pCaloHit->GetInputEnergy();
        accessor[0][2][pixelX][pixelZ] += pCaloHit->GetInputEnergy();
    }

    float minWeight{std::numeric_limits<float>::max()}, maxWeight{-std::numeric_limits<float>::max()};
    for (int r = 0; r < m_imageHeight; ++r)
    {
        for (int c = 0; c < m_imageWidth; ++c)
        {
            // ATTN: Channels are identical, so we only need to look at one
            if (accessor[0][0][r][c] < minWeight)
                minWeight = accessor[0][0][r][c];
            if (accessor[0][0][r][c] > maxWeight)
                maxWeight = accessor[0][0][r][c];
        }
    }
    const float weightRange{maxWeight - minWeight > std::numeric_limits<float>::epsilon() ? maxWeight - minWeight : 1.f};

    int count{0};
    for (int r = 0; r < m_imageHeight; ++r)
    {
        for (int c = 0; c < m_imageWidth; ++c)
        {
            // First we want to get the pixels into the range [0,1] (not [0, 255], transforms.ToTensor handles this)
            accessor[0][0][r][c] = (accessor[0][0][r][c] - minWeight) / weightRange;
            accessor[0][1][r][c] = (accessor[0][1][r][c] - minWeight) / weightRange;
            accessor[0][2][r][c] = (accessor[0][2][r][c] - minWeight) / weightRange;
            if (accessor[0][0][r][c] > 0)
                ++count;

            // Next we want to apply the standard ResNet normalisation (which acts per channel)
            accessor[0][0][r][c] = (accessor[0][0][r][c] - 0.485) / 0.229;
            accessor[0][1][r][c] = (accessor[0][1][r][c] - 0.456) / 0.224;
            accessor[0][2][r][c] = (accessor[0][2][r][c] - 0.406) / 0.225;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoCharacterisationAlgorithm::GetHitRegion(const CaloHitList &caloHitList, double &xMin, double &xMax, double &zMin, double &zMax) const
{
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x(pCaloHit->GetPositionVector().GetX());
        const float z(pCaloHit->GetPositionVector().GetZ());
        if (x < xMin)
            xMin = x;
        if (x > xMax)
            xMax = x;
        if (z < zMin)
            zMin = z;
        if (z > zMax)
            zMax = z;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::ProcessPfoList(const std::string &pfoListName) const
{
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        LArMvaHelper::MvaFeatureVector featureVector;
        CaloHitList caloHitListU, caloHitListV, caloHitListW;
        const HitType viewU{HitType::TPC_VIEW_U}, viewV{HitType::TPC_VIEW_V}, viewW{HitType::TPC_VIEW_W};
        LArPfoHelper::GetCaloHits(pPfo, viewU, caloHitListU);
        LArPfoHelper::GetCaloHits(pPfo, viewV, caloHitListV);
        LArPfoHelper::GetCaloHits(pPfo, viewW, caloHitListW);
        int numGoodViews{0}, numHitsTotal{0};
        for (const CaloHitList &caloHitList : {caloHitListU, caloHitListV, caloHitListW})
        {
            numHitsTotal += caloHitList.size();
            if (caloHitList.size() >= 3)
                ++numGoodViews;
        }
        if ((numGoodViews < 2) || (numHitsTotal < (3 * m_minHitsForGoodView)))
            continue;
        float showerContribution{0.f}, trackContribution{0.f};
        for (const CaloHitList &caloHitList : {caloHitListU, caloHitListV, caloHitListW})
        {
            for (const CaloHit *pCaloHit : caloHitList)
            {
                try
                {
                    const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                    if (!pMCParticle)
                        continue;
                    const int pdg{pMCParticle->GetParticleId()};
                    const int absPdg{std::abs(pdg)};
                    const float energy{pCaloHit->GetInputEnergy()};
                    if (absPdg == E_MINUS || absPdg == PHOTON)
                        showerContribution += energy;
                    else
                        trackContribution += energy;
                }
                catch (const StatusCodeException &)
                {
                }
            }
        }
        const float totalContribution{trackContribution + showerContribution};
        if (totalContribution <= std::numeric_limits<float>::epsilon())
            continue;
        const float trackFraction{trackContribution / totalContribution};
        const int cls{trackFraction > 0.5f ? 1 : 2};
        featureVector.emplace_back(static_cast<double>(cls));
        featureVector.emplace_back(static_cast<double>(trackFraction));
        this->PopulateFeatureVector(caloHitListU, featureVector);
        this->PopulateFeatureVector(caloHitListV, featureVector);
        this->PopulateFeatureVector(caloHitListW, featureVector);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_trainingFileName, true, featureVector));
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListU, "U", RED));
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListV, "V", GREEN));
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListW, "W", BLUE));
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DlPfoCharacterisationAlgorithm::GetTrackFraction(const ParticleFlowObject *const pPfo) const
{
    CaloHitList caloHitListU, caloHitListV, caloHitListW;
    const HitType viewU{HitType::TPC_VIEW_U}, viewV{HitType::TPC_VIEW_V}, viewW{HitType::TPC_VIEW_W};
    LArPfoHelper::GetCaloHits(pPfo, viewU, caloHitListU);
    LArPfoHelper::GetCaloHits(pPfo, viewV, caloHitListV);
    LArPfoHelper::GetCaloHits(pPfo, viewW, caloHitListW);
    int numHitsTotal{0};
    for (const CaloHitList &caloHitList : {caloHitListU, caloHitListV, caloHitListW})
        numHitsTotal += caloHitList.size();
    float showerContribution{0.f}, trackContribution{0.f};
    for (const CaloHitList &caloHitList : {caloHitListU, caloHitListV, caloHitListW})
    {
        for (const CaloHit *pCaloHit : caloHitList)
        {
            try
            {
                const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (!pMCParticle)
                    continue;
                const int pdg{pMCParticle->GetParticleId()};
                const int absPdg{std::abs(pdg)};
                const float energy{pCaloHit->GetInputEnergy()};
                if (absPdg == E_MINUS || absPdg == PHOTON)
                    showerContribution += energy;
                else
                    trackContribution += energy;
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }
    const float totalContribution{trackContribution + showerContribution};
    if (totalContribution <= std::numeric_limits<float>::epsilon())
        return -1.f;

    return trackContribution / totalContribution;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::ChangeCharacterisation(const ParticleFlowObject *const pPfo) const
{
    PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;
    pfoMetadata.m_particleId = pPfo->GetParticleId() == MU_MINUS ? E_MINUS : MU_MINUS;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, pfoMetadata));

    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
    for (const Cluster *const pCluster : clusterList)
    {
        if (pCluster->GetParticleId() == pfoMetadata.m_particleId.Get())
            continue;

        PandoraContentApi::Cluster::Metadata clusterMetadata;
        clusterMetadata.m_particleId = pfoMetadata.m_particleId.Get();
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, clusterMetadata));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoCharacterisationAlgorithm::PopulateFeatureVector(const CaloHitList &caloHitList, LArMvaHelper::MvaFeatureVector &featureVector) const
{
    featureVector.emplace_back(static_cast<double>(caloHitList.size()));
    for (const CaloHit *pCaloHit : caloHitList)
    {
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
        featureVector.emplace_back(static_cast<double>(pCaloHit->GetInputEnergy()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Training", m_training));

    if (m_training)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHitsForGoodView", m_minHitsForGoodView));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_imageHeight));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_imageWidth));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameU", m_modelFileNameU));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameV", m_modelFileNameV));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameW", m_modelFileNameW));

        m_modelFileNameU = LArFileHelper::FindFileInPath(m_modelFileNameU, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileNameU, m_modelU));
        m_modelFileNameV = LArFileHelper::FindFileInPath(m_modelFileNameV, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileNameV, m_modelV));
        m_modelFileNameW = LArFileHelper::FindFileInPath(m_modelFileNameW, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileNameW, m_modelW));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
