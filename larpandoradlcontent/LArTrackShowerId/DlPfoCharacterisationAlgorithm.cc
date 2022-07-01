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
    m_imageHeight{224}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::Run()
{
    if (m_training)
        return this->Train();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::Infer()
{
    PfoList tracksToShowers, showersToTracks;
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({1, 3, m_imageHeight, m_imageWidth}, input);

    const PfoList *pPfoList(nullptr);

    for (const bool inputIsTrack : {true, false})
    {
        const std::string listName{inputIsTrack ? m_trackPfoListName : m_showerPfoListName};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pPfoList));

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            const HitType viewU{HitType::TPC_VIEW_U}, viewV{HitType::TPC_VIEW_V}, viewW{HitType::TPC_VIEW_W};
            CaloHitList caloHitList;
            LArPfoHelper::GetCaloHits(pPfo, viewU, caloHitList);
            LArPfoHelper::GetCaloHits(pPfo, viewV, caloHitList);
            LArPfoHelper::GetCaloHits(pPfo, viewW, caloHitList);
            this->PrepareNetworkInput(caloHitList, input);

            // Run the input through the trained model and get the output accessor
            LArDLHelper::TorchInputVector inputs;
            inputs.emplace_back(input);
            LArDLHelper::TorchOutput output;
            LArDLHelper::Forward(m_model, inputs, output);
            auto outputAccessor = output.accessor<float, 2>();
            // Get class probabilities
            const double output0{exp(outputAccessor[0][0])};
            const double output1{exp(outputAccessor[0][1])};
            const double sum{(output0 + output1) > std::numeric_limits<float>::epsilon() ? output0 + output1 : 1.f};
            const double prob0{output0 / sum}, prob1{output1 / sum};
            if (prob0 > prob1)
            {
                // Shower-like PFO
                if (inputIsTrack)
                {
                    this->ChangeCharacterisation(pPfo);
                    tracksToShowers.emplace_back(pPfo);
                }
            }
            else
            {
                // Track-like PFO
                if (!inputIsTrack)
                {
                    this->ChangeCharacterisation(pPfo);
                    showersToTracks.emplace_back(pPfo);
                }
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

StatusCode DlPfoCharacterisationAlgorithm::Train()
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ProcessPfoList(m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ProcessPfoList(m_showerPfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoCharacterisationAlgorithm::PrepareNetworkInput(const CaloHitList &caloHitList, LArDLHelper::TorchInput &input) const
{
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

    float **weights = new float *[m_imageHeight];
    for (int r = 0; r < m_imageHeight; ++r)
        weights[r] = new float[m_imageWidth]();

    for (const CaloHit *pCaloHit : caloHitList)
    {
        // Determine hit pixel
        const double x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
        const double localX{(x - xMin) / xRange};
        const double localZ{(z - zMin) / zRange};
        const int pixelX{static_cast<int>(std::floor(localX * m_imageWidth))};
        const int pixelZ{static_cast<int>(std::floor(localZ * m_imageHeight))};
        weights[pixelZ][pixelX] += pCaloHit->GetInputEnergy();
    }

    float minWeight{std::numeric_limits<float>::max()}, maxWeight{-std::numeric_limits<float>::max()};
    for (int r = 0; r < m_imageHeight; ++r)
    {
        for (int c = 0; c < m_imageWidth; ++c)
        {
            if (weights[r][c] < minWeight)
                minWeight = weights[r][c];
            if (weights[r][c] > maxWeight)
                maxWeight = weights[r][c];
        }
    }
    const float weightRange{maxWeight - minWeight > std::numeric_limits<float>::epsilon() ? maxWeight - minWeight : 1.f};

    auto accessor = input.accessor<float, 4>();
    // Channel mapping is UVW => RGB (012)
    const int channel{caloHitList.front()->GetHitType() - HitType::TPC_VIEW_U};
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const double x(pCaloHit->GetPositionVector().GetX());
        const double z(pCaloHit->GetPositionVector().GetZ());
        const double localX{(x - xMin) / xRange};
        const double localZ{(z - zMin) / zRange};
        const int pixelX{static_cast<int>(std::floor(localX * m_imageWidth))};
        const int pixelZ{static_cast<int>(std::floor(localZ * m_imageHeight))};
        accessor[0][channel][pixelZ][pixelX] = (weights[pixelZ][pixelX] - minWeight) / weightRange;
    }

    for (int r = 0; r < m_imageHeight; ++r)
        delete [] weights[r];
    delete[] weights;
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
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        LArMvaHelper::MvaFeatureVector featureVector;
        float trackContribution{0.f}, showerContribution{0.f};
        CaloHitList caloHitListU, caloHitListV, caloHitListW;
        const HitType viewU{HitType::TPC_VIEW_U}, viewV{HitType::TPC_VIEW_V}, viewW{HitType::TPC_VIEW_W};
        LArPfoHelper::GetCaloHits(pPfo, viewU, caloHitListU);
        LArPfoHelper::GetCaloHits(pPfo, viewV, caloHitListV);
        LArPfoHelper::GetCaloHits(pPfo, viewW, caloHitListW);
        int numGoodViews{0};
        if (caloHitListU.size() >= 5)
            ++numGoodViews;
        if (caloHitListV.size() >= 5)
            ++numGoodViews;
        if (caloHitListW.size() >= 5)
            ++numGoodViews;
        if (numGoodViews < 2 || (caloHitListU.size() + caloHitListV.size() + caloHitListW.size() < 15))
            continue;
        for (const CaloHitList &caloHitList : {caloHitListU, caloHitListV, caloHitListW})
        {
            for (const CaloHit *pCaloHit : caloHitList)
            {
                try
                {
                    const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                    const int pdg{std::abs(pMCParticle->GetParticleId())};
                    if (pdg == E_MINUS || pdg == PHOTON)
                        showerContribution += pCaloHit->GetInputEnergy();
                    else
                        trackContribution += pCaloHit->GetInputEnergy();
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
        featureVector.emplace_back(static_cast<double>(trackFraction));
        this->PopulateFeatureVector(caloHitListU, featureVector);
        this->PopulateFeatureVector(caloHitListV, featureVector);
        this->PopulateFeatureVector(caloHitListW, featureVector);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_trainingFileName,
            trackFraction > 0.5f, featureVector));
    }

    return STATUS_CODE_SUCCESS;
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
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_imageHeight));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_imageWidth));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", m_modelFileName));
        m_modelFileName = LArFileHelper::FindFileInPath(m_modelFileName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileName, m_model));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
