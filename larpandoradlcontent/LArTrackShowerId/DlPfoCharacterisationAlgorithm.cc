/**
 *  @file   larpandoracontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArCalorimetryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/QuadTree.h"

#include "larpandoradlcontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.h"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlPfoCharacterisationAlgorithm::DlPfoCharacterisationAlgorithm() :
    m_visualize{false},
    m_training{false}
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
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }

    PfoList tracksToShowers, showersToTracks;

    const PfoList *pPfoList(nullptr);
    for (const bool inputIsTrack : {true, false})
    {
        const std::string listName{inputIsTrack ? m_trackPfoListName : m_showerPfoListName};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pPfoList));

        if (pPfoList->empty())
            continue;

        LArDLHelper::TorchInput input;
        LArDLHelper::InitialiseInput({static_cast<int>(pPfoList->size()), 243}, input);
        this->PrepareNetworkInput(*pPfoList, input);

        // Run the input through the trained model and get the output accessor
        LArDLHelper::TorchInputVector inputs;
        inputs.emplace_back(input);
        LArDLHelper::TorchOutput output;
        LArDLHelper::Forward(m_model, inputs, output);
        auto outputAccessor = output.accessor<float, 2>();
        int i{0};
        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            // Get class probabilities
            const double outputTrack{exp(outputAccessor[i][0])};
            const double outputShower{exp(outputAccessor[i][1])};
            const double sum{(outputShower + outputTrack) > std::numeric_limits<float>::epsilon() ? outputShower + outputTrack : 1.f};
            const double probShower{outputShower / sum}, probTrack{outputTrack / sum};

            if (inputIsTrack && probShower > probTrack)
            {
                this->ChangeCharacterisation(pPfo);
                tracksToShowers.emplace_back(pPfo);
            }
            else if (!inputIsTrack && probTrack > probShower)
            {
                this->ChangeCharacterisation(pPfo);
                showersToTracks.emplace_back(pPfo);
            }

            if (m_visualize)
            {
                CaloHitList caloHitListU, caloHitListV, caloHitListW;
                LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_U, caloHitListU);
                LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_V, caloHitListV);
                LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_W, caloHitListW);
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListU, "U", probTrack > probShower ? BLUE : RED));
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListV, "V", probTrack > probShower ? BLUE : RED));
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListW, "W", probTrack > probShower ? BLUE : RED));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
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
    this->ProcessPfoList(m_trackPfoListName);
    this->ProcessPfoList(m_showerPfoListName);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoCharacterisationAlgorithm::PrepareNetworkInput(const PfoList& pfoList, LArDLHelper::TorchInput &input) const
{
    int p{0};
    auto accessor = input.accessor<float, 2>();
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        int e{0};
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            CaloHitList caloHitList;
            LArPfoHelper::GetCaloHits(pPfo, view, caloHitList);

            FloatVector longitudinalProfile, transverseProfile;
            try
            {
                CartesianVector origin(0.f, 0.f, 0.f), longitudinal(0.f, 0.f, 0.f), transverse(0.f, 0.f, 0.f);
                LArCalorimetryHelper::GetPrincipalAxes(caloHitList, origin, longitudinal, transverse);
                longitudinalProfile = LArCalorimetryHelper::GetLongitudinalAdcProfile(caloHitList, origin, longitudinal, 0.1f, 14.f);
                transverseProfile = LArCalorimetryHelper::GetTransverseAdcProfile(caloHitList, origin, transverse, 0.15f, 9.f);
            }
            catch (...)
            {
            }

            FloatVector startProfile(30), endProfile(30);
            int i{0};
            // Starting profile
            for (const float adc : longitudinalProfile)
            {
                startProfile[i] = adc;
                ++i;
                if (i == 30)
                    break;
            }
            // Ending profile
            i = 29;
            for (auto iter = longitudinalProfile.rbegin(); iter != longitudinalProfile.rend(); ++iter)
            {
                endProfile[i] = *iter;
                --i;
                if (i < 0)
                    break;
            }
            // Transverse profile
            size_t length{transverseProfile.size()};
            size_t midPoint{length / 2};
            const int nRadialBins{10};
            FloatVector tProfile(2 * nRadialBins + 1);
            for (i = -nRadialBins; i <= nRadialBins; ++i)
            {
                if ((i + static_cast<long>(midPoint)) < 0)
                    continue;
                if ((i + static_cast<long>(midPoint)) >= static_cast<long>(transverseProfile.size()))
                    break;
                tProfile[i + nRadialBins] = transverseProfile[i + midPoint];
            }

            for (float deposit : startProfile)
                accessor[p][e++] = deposit;
            for (float deposit : endProfile)
                accessor[p][e++] = deposit;
            for (float deposit : tProfile)
                accessor[p][e++] = deposit;
        }
        ++p;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::ProcessPfoList(const std::string &pfoListName) const
{
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));
    this->MakeTrainingImage(*pPfoList);

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        CaloHitList caloHitList, caloHitListU, caloHitListV, caloHitListW;
        LArPfoHelper::GetAllCaloHits(pPfo, caloHitList);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_U, caloHitListU);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_V, caloHitListV);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_W, caloHitListW);
        int nGoodViews{0};
        nGoodViews += caloHitListU.size() >= 15 ? 1 : 0;
        nGoodViews += caloHitListV.size() >= 15 ? 1 : 0;
        nGoodViews += caloHitListW.size() >= 15 ? 1 : 0;

        if (caloHitList.size() < 15 || nGoodViews < 2)
            continue;

        float showerContribution{0.f}, trackContribution{0.f}, neutronContribution{0.f}, totalContribution{0.f};
        for (const CaloHit *pCaloHit : caloHitList)
        {
            try
            {
                const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (!pMCParticle)
                    continue;
                const int absPdg{static_cast<int>(std::abs(pMCParticle->GetParticleId()))};
                const float adc{pCaloHit->GetInputEnergy()};
                const MCParticle *pParent{pMCParticle};
                while (!pParent->GetParentList().empty())
                {
                    pParent = pParent->GetParentList().front();
                    if (std::abs(pParent->GetParticleId()) == NEUTRON)
                    {
                        neutronContribution += adc;
                        break;
                    }
                }
                totalContribution += adc;
                if (absPdg == E_MINUS || absPdg == PHOTON)
                    showerContribution += adc;
                else
                    trackContribution += adc;
            }
            catch (const StatusCodeException &)
            {
            }
        }
        if (totalContribution <= std::numeric_limits<float>::epsilon())
            continue;
        const float trackFraction{trackContribution / totalContribution};
        const float neutronFraction{neutronContribution / totalContribution};
        const bool isDownstreamNeutron{neutronFraction > 0.5f ? true : false};
        if (isDownstreamNeutron)
            continue;

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(trackFraction));
        const int cls{trackFraction > 0.5f ? 1 : 2};
        featureVector.emplace_back(static_cast<double>(cls));

        this->ProcessPfoView(caloHitListU, featureVector);
        this->ProcessPfoView(caloHitListV, featureVector);
        this->ProcessPfoView(caloHitListW, featureVector);
        // Temporarily stop producing calorimetric output. Need to figure out how to combine this with image info
        // Probably need to alter image production to process one PFO-at-a-time and pass in global info from a preprocessing
        // step to ensure everything can be added to the same feature vector
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_trainingFileName, true, featureVector));

        //std::cout << "NHits " << caloHitList.size() << " Neutron " << isDownstreamNeutron << " Track Frac " << trackFraction << " CLS " << cls << std::endl;
        //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "PFO", BLUE));
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlPfoCharacterisationAlgorithm::ProcessPfoView(const CaloHitList &caloHitList, LArMvaHelper::MvaFeatureVector &featureVector) const
{
    CartesianVector origin(0.f, 0.f, 0.f), longitudinal(0.f, 0.f, 0.f), transverse(0.f, 0.f, 0.f);
    try
    {
        LArCalorimetryHelper::GetPrincipalAxes(caloHitList, origin, longitudinal, transverse);
        FloatVector longitudinalProfile{LArCalorimetryHelper::GetLongitudinalAdcProfile(caloHitList, origin, longitudinal, 0.1f, 14.f)};
        FloatVector transverseProfile{LArCalorimetryHelper::GetTransverseAdcProfile(caloHitList, origin, transverse, 0.15f, 9.f)};
        this->PopulateFeatureVector(longitudinalProfile, transverseProfile, featureVector);
    }
    catch (...)
    {
        // Still need to produce empty profiles if a view has insuffient hits for PCA
        FloatVector longitudinalProfile, transverseProfile;
        this->PopulateFeatureVector(longitudinalProfile, transverseProfile, featureVector);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::MakeTrainingImage(const PfoList &pfoList) const
{
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }

    CaloHitList caloHitListU, caloHitListV, caloHitListW;
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_U, caloHitListU);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_V, caloHitListV);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_W, caloHitListW);
    }
    QuadTree quadTreeU(caloHitListU);
    QuadTree quadTreeV(caloHitListV);
    QuadTree quadTreeW(caloHitListW);

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        LArMvaHelper::MvaFeatureVector featureVector;
        CaloHitList pfoHitListU, pfoHitListV, pfoHitListW, pfoHitList3D, pfoHitListAll;
        LArPfoHelper::GetAllCaloHits(pPfo, pfoHitListAll);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_U, pfoHitListU);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_V, pfoHitListV);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_VIEW_W, pfoHitListW);
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_3D, pfoHitList3D);
        if (pfoHitList3D.empty())
            continue;

        int nGoodViews{0};
        nGoodViews += pfoHitListU.size() >= 15 ? 1 : 0;
        nGoodViews += pfoHitListV.size() >= 15 ? 1 : 0;
        nGoodViews += pfoHitListW.size() >= 15 ? 1 : 0;

        if (pfoHitListAll.size() < 15 || nGoodViews < 2)
            continue;

        float showerContribution{0.f}, trackContribution{0.f}, neutronContribution{0.f}, totalContribution{0.f};
        for (const CaloHit *pCaloHit : pfoHitListAll)
        {
            try
            {
                const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                if (!pMCParticle)
                    continue;
                const int absPdg{static_cast<int>(std::abs(pMCParticle->GetParticleId()))};
                const float adc{pCaloHit->GetInputEnergy()};
                const MCParticle *pParent{pMCParticle};
                while (!pParent->GetParentList().empty())
                {
                    pParent = pParent->GetParentList().front();
                    if (std::abs(pParent->GetParticleId()) == NEUTRON)
                    {
                        neutronContribution += adc;
                        break;
                    }
                }
                totalContribution += adc;
                if (absPdg == E_MINUS || absPdg == PHOTON)
                    showerContribution += adc;
                else
                    trackContribution += adc;
            }
            catch (const StatusCodeException &)
            {
            }
        }
        if (totalContribution <= std::numeric_limits<float>::epsilon())
            continue;
        const float trackFraction{trackContribution / totalContribution};
        const float neutronFraction{neutronContribution / totalContribution};
        const bool isDownstreamNeutron{neutronFraction > 0.5f ? true : false};
        if (isDownstreamNeutron)
            continue;

        featureVector.emplace_back(static_cast<double>(trackFraction));
        const int cls{trackFraction > 0.5f ? 1 : 2};
        featureVector.emplace_back(static_cast<double>(cls));

        float xMin{std::numeric_limits<float>::max()}, xMax{-std::numeric_limits<float>::max()};
        float zMinU{std::numeric_limits<float>::max()}, zMaxU{-std::numeric_limits<float>::max()};
        for (const CaloHit *pCaloHit : pfoHitListU)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            xMin = std::min(pos.GetX(), xMin);
            xMax = std::max(pos.GetX(), xMax);
            zMinU = std::min(pos.GetZ(), zMinU);
            zMaxU = std::max(pos.GetZ(), zMaxU);
        }
        float zMinV{std::numeric_limits<float>::max()}, zMaxV{-std::numeric_limits<float>::max()};
        for (const CaloHit *pCaloHit : pfoHitListV)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            xMin = std::min(pos.GetX(), xMin);
            xMax = std::max(pos.GetX(), xMax);
            zMinV = std::min(pos.GetZ(), zMinV);
            zMaxV = std::max(pos.GetZ(), zMaxV);
        }
        float zMinW{std::numeric_limits<float>::max()}, zMaxW{-std::numeric_limits<float>::max()};
        for (const CaloHit *pCaloHit : pfoHitListW)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            xMin = std::min(pos.GetX(), xMin);
            xMax = std::max(pos.GetX(), xMax);
            zMinW = std::min(pos.GetZ(), zMinW);
            zMaxW = std::max(pos.GetZ(), zMaxW);
        }
        CartesianPointVector pointVector;
        for (const CaloHit *pCaloHit : pfoHitList3D)
            pointVector.emplace_back(pCaloHit->GetPositionVector());

        ClusterList clusters3D;
        LArPfoHelper::GetClusters(pPfo, TPC_3D, clusters3D);
        CartesianVector vertex(0, 0, 0), endpoint(0, 0, 0);
        const Cluster* pCluster{clusters3D.front()};
        LArClusterHelper::GetExtremalCoordinates(pCluster, vertex, endpoint);
        const LArTPC *const pLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
        const float wirePitch(pLArTPC->GetWirePitchW());
        LArTrackStateVector trackTraj;
        LArPfoHelper::GetSlidingFitTrajectory(pointVector, vertex, 20.f, wirePitch, trackTraj);
        if (!trackTraj.empty())
        {
            const CartesianVector &direction{(trackTraj.begin())->GetDirection()};
            const CartesianVector &dirU{LArGeometryHelper::ProjectDirection(this->GetPandora(), direction, TPC_VIEW_U)};
            const CartesianVector &dirV{LArGeometryHelper::ProjectDirection(this->GetPandora(), direction, TPC_VIEW_V)};
            const CartesianVector &dirW{LArGeometryHelper::ProjectDirection(this->GetPandora(), direction, TPC_VIEW_W)};

            // We want a 64 x 64 region, so if it's too small expand it in the downstream direction and if it's too large maintain the upstream portion
            if (direction.GetX() < 0.f)
                xMin = xMax - 64.f;
            else
                xMax = xMin + 64.f;
            xMin -= 0.5f; xMax += 0.5f;
            if (!pfoHitListU.empty())
            {
                if (dirU.GetZ() < 0.f)
                    zMinU = zMaxU - 64.f;
                else
                    zMaxU = zMinU + 64.f;
                zMinU -= 0.5f; zMaxU += 0.5f;
            }
            if (!pfoHitListV.empty())
            {
                if (dirV.GetZ() < 0.f)
                    zMinV = zMaxV - 64.f;
                else
                    zMaxV = zMinV + 64.f;
                zMinV -= 0.5f; zMaxV += 0.5f;
            }
            if (!pfoHitListW.empty())
            {
                if (dirW.GetZ() < 0.f)
                    zMinW = zMaxW - 64.f;
                else
                    zMaxW = zMinW + 64.f;
                zMinW -= 0.5f; zMaxW += 0.5f;
            }

            CaloHitList backgroundHitListU;
            if (!pfoHitListU.empty())
                quadTreeU.Find(xMin, zMinU, xMax, zMaxU, backgroundHitListU);

            CaloHitList backgroundHitListV;
            if (!pfoHitListV.empty())
                quadTreeV.Find(xMin, zMinV, xMax, zMaxV, backgroundHitListV);
           
            CaloHitList backgroundHitListW;
                quadTreeW.Find(xMin, zMinW, xMax, zMaxW, backgroundHitListW);

            CaloHitList retainedPfoHitListU, retainedPfoHitListV, retainedPfoHitListW;
            CaloHitList rejectedPfoHitListU, rejectedPfoHitListV, rejectedPfoHitListW;
            for (const CaloHit *const pCaloHit : pfoHitListU)
            {
                const CartesianVector &pos{pCaloHit->GetPositionVector()};
                if (xMin <= pos.GetX() && pos.GetX() <= xMax && zMinU <= pos.GetZ() && pos.GetZ() <= zMaxU)
                    retainedPfoHitListU.emplace_back(pCaloHit);
                else
                    rejectedPfoHitListU.emplace_back(pCaloHit);
            }
            for (const CaloHit *const pCaloHit : pfoHitListV)
            {
                const CartesianVector &pos{pCaloHit->GetPositionVector()};
                if (xMin <= pos.GetX() && pos.GetX() <= xMax && zMinV <= pos.GetZ() && pos.GetZ() <= zMaxV)
                    retainedPfoHitListV.emplace_back(pCaloHit);
                else
                    rejectedPfoHitListV.emplace_back(pCaloHit);
            }
            for (const CaloHit *const pCaloHit : pfoHitListW)
            {
                const CartesianVector &pos{pCaloHit->GetPositionVector()};
                if (xMin <= pos.GetX() && pos.GetX() <= xMax && zMinW <= pos.GetZ() && pos.GetZ() <= zMaxW)
                    retainedPfoHitListW.emplace_back(pCaloHit);
                else
                    rejectedPfoHitListW.emplace_back(pCaloHit);
            }

            for (auto iter = backgroundHitListU.begin(); iter != backgroundHitListU.end(); )
            {
                if (std::find(retainedPfoHitListU.begin(), retainedPfoHitListU.end(), *iter) != retainedPfoHitListU.end())
                    iter = backgroundHitListU.erase(iter);
                else
                    ++iter;
            }
            for (auto iter = backgroundHitListV.begin(); iter != backgroundHitListV.end(); )
            {
                if (std::find(retainedPfoHitListV.begin(), retainedPfoHitListV.end(), *iter) != retainedPfoHitListV.end())
                    iter = backgroundHitListV.erase(iter);
                else
                    ++iter;
            }
            for (auto iter = backgroundHitListW.begin(); iter != backgroundHitListW.end(); )
            {
                if (std::find(retainedPfoHitListW.begin(), retainedPfoHitListW.end(), *iter) != retainedPfoHitListW.end())
                    iter = backgroundHitListW.erase(iter);
                else
                    ++iter;
            }

            featureVector.emplace_back(static_cast<double>(xMin));
            featureVector.emplace_back(static_cast<double>(xMax));
            featureVector.emplace_back(static_cast<double>(zMinU));
            featureVector.emplace_back(static_cast<double>(zMaxU));
            featureVector.emplace_back(static_cast<double>(zMinV));
            featureVector.emplace_back(static_cast<double>(zMaxV));
            featureVector.emplace_back(static_cast<double>(zMinW));
            featureVector.emplace_back(static_cast<double>(zMaxW));
            for (const CaloHitList &retainedPfoHitList : {retainedPfoHitListU, retainedPfoHitListV, retainedPfoHitListW})
            {
                featureVector.emplace_back(static_cast<double>(retainedPfoHitList.size()));
                for (const CaloHit *const pCaloHit : retainedPfoHitList)
                {
                    const CartesianVector &pos{pCaloHit->GetPositionVector()};
                    featureVector.emplace_back(static_cast<double>(pos.GetX()));
                    featureVector.emplace_back(static_cast<double>(pos.GetZ()));
                    featureVector.emplace_back(static_cast<double>(pCaloHit->GetInputEnergy()));
                }
            }
            for (const CaloHitList &backgroundHitList : {backgroundHitListU, backgroundHitListV, backgroundHitListW})
            {
                featureVector.emplace_back(static_cast<double>(backgroundHitList.size()));
                for (const CaloHit *const pCaloHit : backgroundHitList)
                {
                    const CartesianVector &pos{pCaloHit->GetPositionVector()};
                    featureVector.emplace_back(static_cast<double>(pos.GetX()));
                    featureVector.emplace_back(static_cast<double>(pos.GetZ()));
                    featureVector.emplace_back(static_cast<double>(-pCaloHit->GetInputEnergy()));
                }
            }
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_trainingFileName, true, featureVector));

            if (m_visualize)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &retainedPfoHitListW, "PFO", BLUE));
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &rejectedPfoHitListW, "Rejects", GREEN));
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &backgroundHitListW, "Background", RED));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }
        }

        //const LArShowerPCA &showerPca{LArPfoHelper::GetPrincipalComponents(pPfo, pVertex)};
        //const CartesianVector &showerDirection{showerPca.GetPrimaryAxis()};
        //std::cout << "Shower Direction " << showerDirection << std::endl;
        // Plot the PFO in one colour, surrounding hits in another and the full PFO ina third so we can see what decisions are being made
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

void DlPfoCharacterisationAlgorithm::PopulateFeatureVector(const FloatVector &longitudinalProfile, const FloatVector &transverseProfile,
    LArMvaHelper::MvaFeatureVector &featureVector) const
{
    FloatVector startProfile(30), endProfile(30);
    int i{0};
    // Starting profile
    for (const float adc : longitudinalProfile)
    {
        startProfile[i] = adc;
        ++i;
        if (i == 30)
            break;
    }
    // Ending profile
    i = 29;
    for (auto iter = longitudinalProfile.rbegin(); iter != longitudinalProfile.rend(); ++iter)
    {
        endProfile[i] = *iter;
        --i;
        if (i < 0)
            break;
    }
    // Transverse profile
    size_t length{transverseProfile.size()};
    size_t midPoint{length / 2};
    const int nRadialBins{10};
    FloatVector tProfile(2 * nRadialBins + 1);
    for (i = -nRadialBins; i <= nRadialBins; ++i)
    {
        if ((i + static_cast<long>(midPoint)) < 0)
            continue;
        if ((i + static_cast<long>(midPoint)) >= static_cast<long>(transverseProfile.size()))
            break;
        tProfile[i + nRadialBins] = transverseProfile[i + midPoint];
    }

    featureVector.emplace_back(static_cast<double>(startProfile.size()));
    for (const float adc : startProfile)
        featureVector.emplace_back(static_cast<double>(adc));
    featureVector.emplace_back(static_cast<double>(endProfile.size()));
    for (const float adc : endProfile)
        featureVector.emplace_back(static_cast<double>(adc));
    featureVector.emplace_back(static_cast<double>(tProfile.size()));
    for (const float adc : tProfile)
        featureVector.emplace_back(static_cast<double>(adc));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Training", m_training));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    if (m_training)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", m_modelFileName));

        m_modelFileName = LArFileHelper::FindFileInPath(m_modelFileName, "FW_SEARCH_PATH");
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_modelFileName, m_model));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
