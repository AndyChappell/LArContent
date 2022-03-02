/**
 *  @file   larpandoracontent/LArTrackShowerId/HitTruthTaggingAlgorithm.cc
 *
 *  @brief  Implementation of the branch growing algorithm base class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTrackShowerId/HitTruthTaggingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

enum HitTruthTaggingAlgorithm::Tag : int { SHOWER = 1, TRACK };

HitTruthTaggingAlgorithm::HitTruthTaggingAlgorithm() : m_minTrackRatio{0.5f}, m_visualize{false}, m_writeFeatures{true},
    m_confusionU(), m_confusionV(), m_confusionW()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitTruthTaggingAlgorithm::~HitTruthTaggingAlgorithm()
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_true_shower", m_confusionU[1][1]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_false_shower", m_confusionU[0][1]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_false_track", m_confusionU[1][0]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "u_true_track", m_confusionU[0][0]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_true_shower", m_confusionV[1][1]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_false_shower", m_confusionV[0][1]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_false_track", m_confusionV[1][0]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "v_true_track", m_confusionV[0][0]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_true_shower", m_confusionW[1][1]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_false_shower", m_confusionW[0][1]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_false_track", m_confusionW[1][0]));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "confusion_tree", "w_true_track", m_confusionW[0][0]));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "confusion_tree"));
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "confusion_tree", "confusion.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "HitTruthTaggingAlgorithm: Unable to write confusion_tree to file" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitTruthTaggingAlgorithm::TagMuonClusters(CaloHitList &trackHitList, CaloHitList &showerHitList) const
{
    for (std::string clusterListName : m_muonClusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        if (PandoraContentApi::GetList(*this, clusterListName, pClusterList) != STATUS_CODE_SUCCESS)
            continue;

        std::map<const MCParticle *, size_t> mcToLargestClusterMap;
        std::map<const MCParticle *, const LArPointingCluster *> mcToPointingClusterMap;
        for (const Cluster *const pCluster : *pClusterList)
        {
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
            caloHitList.insert(caloHitList.end(), isolatedHits.begin(), isolatedHits.end());
            if (!caloHitList.empty())
            {


                const MCParticle *pMCParticle{nullptr};

                bool foundMC{false};
                for (const CaloHit *const pRefCaloHit : caloHitList)
                {
                    try
                    {
                        pMCParticle = MCParticleHelper::GetMainMCParticle(pRefCaloHit);
                        foundMC = true;
                        break;
                    }
                    catch (const StatusCodeException &)
                    {
                    }
                }
                if (!foundMC)
                    continue;

                if (mcToLargestClusterMap.find(pMCParticle) != mcToLargestClusterMap.end())
                {
                    if (caloHitList.size() > mcToLargestClusterMap[pMCParticle])
                    {
                        mcToLargestClusterMap[pMCParticle] = caloHitList.size();
                        if (caloHitList.size() >= 5)
                        {
                            if (mcToPointingClusterMap.find(pMCParticle) != mcToPointingClusterMap.end())
                            {
                                delete mcToPointingClusterMap[pMCParticle];
                            }
                            mcToPointingClusterMap[pMCParticle] = new LArPointingCluster(pCluster);
                        }
                    }
                }
                else
                {
                    mcToLargestClusterMap[pMCParticle] = caloHitList.size();
                    if (caloHitList.size() >= 5)
                    {
                        mcToPointingClusterMap[pMCParticle] = new LArPointingCluster(pCluster);
                    }
                }
            }
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
            caloHitList.insert(caloHitList.end(), isolatedHits.begin(), isolatedHits.end());

            if (!caloHitList.empty())
            {
                int tag{TRACK};
                const MCParticle *pMCParticle{nullptr};

                bool foundMC{false};
                for (const CaloHit *const pRefCaloHit : caloHitList)
                {
                    try
                    {
                        pMCParticle = MCParticleHelper::GetMainMCParticle(pRefCaloHit);
                        foundMC = true;
                        break;
                    }
                    catch (const StatusCodeException &)
                    {
                    }
                }
                if (!foundMC)
                    continue;
                const LArMCParticle *const pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
                if (!pLArMCParticle)
                    continue;

                if (LArMCParticleHelper::IsCapture(pLArMCParticle) || LArMCParticleHelper::IsNuclear(pLArMCParticle) ||
                    LArMCParticleHelper::IsIonisation(pLArMCParticle))
                {
                    tag = SHOWER;
                }
                else if(LArMCParticleHelper::IsInelasticScatter(pLArMCParticle))
                {
                    if (caloHitList.size() >= std::max(5.f, m_minTrackRatio * mcToLargestClusterMap[pMCParticle]))
                        tag = TRACK;
                    else
                        tag = SHOWER;
                }
                else
                {
                    if (caloHitList.size() >= 5)
                    {
                        const LArPointingCluster thisPointingCluster(pCluster);
                        const LArPointingCluster &mainPointingCluster{*mcToPointingClusterMap[pMCParticle]};
                        if (thisPointingCluster.GetCluster() == mainPointingCluster.GetCluster())
                        {
                            tag = TRACK;
                        }
                        else
                        {
                            const LArPointingCluster::Vertex *pTheseVertices[]{&thisPointingCluster.GetInnerVertex(),
                                &thisPointingCluster.GetOuterVertex()};
                            const LArPointingCluster::Vertex *pMainVertices[]{&mainPointingCluster.GetInnerVertex(),
                                &mainPointingCluster.GetOuterVertex()};
                            float minDistance{std::numeric_limits<float>::max()};
                            const LArPointingCluster::Vertex *pMainVertex{nullptr};

                            for (const LArPointingCluster::Vertex *pVertex1 : pTheseVertices)
                            {
                                for (const LArPointingCluster::Vertex *pVertex2 : pMainVertices)
                                {
                                    const float distance{(pVertex1->GetPosition() - pVertex2->GetPosition()).GetMagnitudeSquared()};
                                    if (distance < minDistance) 
                                    {
                                        minDistance = distance;
                                        pMainVertex = pVertex2;
                                    }
                                }
                            }

                            const CartesianVector &mainDirection{pMainVertex->GetDirection()};
                            CaloHitList mainHitList;
                            mainPointingCluster.GetCluster()->GetOrderedCaloHitList().FillCaloHitList(mainHitList);
                            FloatVector mainProj, clusterProj;
                            const CartesianVector &vertex{pMainVertex->GetPosition()};
                            for (const CaloHit *pHit : mainHitList)
                            {
                                const CartesianVector localDirection{pHit->GetPositionVector() - vertex};
                                mainProj.emplace_back(localDirection.GetDotProduct(mainDirection));
                            }
                            for (const CaloHit *pHit : caloHitList)
                            {
                                const CartesianVector localDirection{pHit->GetPositionVector() - vertex};
                                clusterProj.emplace_back(localDirection.GetDotProduct(mainDirection));
                            }
                            const float clusterMin{*std::min_element(clusterProj.begin(), clusterProj.end())},
                                clusterMax{*std::max_element(clusterProj.begin(), clusterProj.end())};
                            const float mainMin{*std::min_element(mainProj.begin(), mainProj.end())},
                                mainMax{*std::max_element(mainProj.begin(), mainProj.end())};
                            if ((clusterMin > (mainMin + .5f) && clusterMin < (mainMax - 0.5f)) ||
                                (clusterMax > (mainMin + .5f) && clusterMax < (mainMax - 0.5f)) ||
                                (clusterMin < mainMin && clusterMax > mainMax))
                            {
                                tag = SHOWER;
                            }
                            else
                            {
                                tag = TRACK;
                            }
                        }
                    }
                    else
                    {
                        tag = SHOWER;
                    }
                }

                for (const CaloHit *const pHit : caloHitList)
                {
                    switch (tag)
                    {
                        case TRACK:
                            trackHitList.emplace_back(pHit);
                            break;
                        case SHOWER:
                            showerHitList.emplace_back(pHit);
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitTruthTaggingAlgorithm::TagNonMuonClusters(CaloHitList &trackHitList, CaloHitList &showerHitList) const
{
    for (std::string clusterListName : m_nonMuonClusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        if (PandoraContentApi::GetList(*this, clusterListName, pClusterList) != STATUS_CODE_SUCCESS)
            continue;

        for (const Cluster *const pCluster : *pClusterList)
        {
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
            caloHitList.insert(caloHitList.end(), isolatedHits.begin(), isolatedHits.end());

            if (!caloHitList.empty())
            {
                try
                {
                    int tag{TRACK};
                    const CaloHit *pCaloHit{caloHitList.front()};
                    const MCParticle *const pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                    const LArMCParticle *const pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
                    if (!pLArMCParticle)
                        continue;
                    const int pdg{std::abs(pMCParticle->GetParticleId())};

                    if (LArMCParticleHelper::IsCapture(pLArMCParticle) || LArMCParticleHelper::IsNuclear(pLArMCParticle) ||
                        LArMCParticleHelper::IsIonisation(pLArMCParticle))
                    {
                        tag = SHOWER;
                    }
                    else if(LArMCParticleHelper::IsInelasticScatter(pLArMCParticle))
                    {
                        if (caloHitList.size() >= 5)
                        {
                            if (pdg == PHOTON || pdg == E_MINUS)
                                tag = SHOWER;
                            else
                                tag = TRACK;
                        }
                        else
                        {
                            tag = SHOWER;
                        }
                    }
                    else if (pdg == PHOTON)
                    {
                        tag = SHOWER;
                    }
                    else if (pdg == E_MINUS)
                    {
                        tag = SHOWER;
                    }
                    else
                    {
                        tag = TRACK;
                    }
                    for (const CaloHit *const pHit : caloHitList)
                    {
                        switch (tag)
                        {
                            case TRACK:
                                trackHitList.emplace_back(pHit);
                                break;
                            case SHOWER:
                                showerHitList.emplace_back(pHit);
                                break;
                            default:
                                break;
                        }
                    }
                }
                catch (const StatusCodeException &)
                {
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitTruthTaggingAlgorithm::OutputTruthTagging(const CaloHitList &trackHitList, const CaloHitList &showerHitList) const
{
    LArMvaHelper::MvaFeatureVector featureVectorU, featureVectorV, featureVectorW;

    for (const CaloHit *pCaloHit : trackHitList)
    {
        LArMvaHelper::MvaFeatureVector *pFeatureVector{nullptr};
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                pFeatureVector = &featureVectorU;
                break;
            case TPC_VIEW_V:
                pFeatureVector = &featureVectorV;
                break;
            case TPC_VIEW_W:
                pFeatureVector = &featureVectorW;
                break;
            default:
                continue;
        }

        pFeatureVector->emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
        pFeatureVector->emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
        pFeatureVector->emplace_back(static_cast<double>(TRACK));
        pFeatureVector->emplace_back(static_cast<double>(pCaloHit->GetInputEnergy()));
    }

    for (const CaloHit *pCaloHit : showerHitList)
    {
        LArMvaHelper::MvaFeatureVector *pFeatureVector{nullptr};
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                pFeatureVector = &featureVectorU;
                break;
            case TPC_VIEW_V:
                pFeatureVector = &featureVectorV;
                break;
            case TPC_VIEW_W:
                pFeatureVector = &featureVectorW;
                break;
            default:
                continue;
        }

        pFeatureVector->emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
        pFeatureVector->emplace_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
        pFeatureVector->emplace_back(static_cast<double>(SHOWER));
        pFeatureVector->emplace_back(static_cast<double>(pCaloHit->GetInputEnergy()));
    }

    // Add number of hits to end of vector than rotate (more efficient than direct insert at front)
    featureVectorU.push_back(static_cast<double>(featureVectorU.size() / 4));
    std::rotate(featureVectorU.rbegin(), featureVectorU.rbegin() + 1, featureVectorU.rend());
    featureVectorV.push_back(static_cast<double>(featureVectorV.size() / 4));
    std::rotate(featureVectorV.rbegin(), featureVectorV.rbegin() + 1, featureVectorV.rend());
    featureVectorW.push_back(static_cast<double>(featureVectorW.size() / 4));
    std::rotate(featureVectorW.rbegin(), featureVectorW.rbegin() + 1, featureVectorW.rend());

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_outputFileName + "U.csv", true, featureVectorU));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_outputFileName + "V.csv", true, featureVectorV));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArMvaHelper::ProduceTrainingExample(m_outputFileName + "W.csv", true, featureVectorW));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitTruthTaggingAlgorithm::ValidateInference(const CaloHitList &trackHitList, const CaloHitList &showerHitList)
{
    if (m_visualize)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }

    CaloHitList correctHitListW, incorrectHitListW;
    for (const CaloHit *pCaloHit : trackHitList)
    {
        HitType view{pCaloHit->GetHitType()};
        const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        const float pTrack{pLArCaloHit->GetTrackProbability()};
        const float pShower{pLArCaloHit->GetShowerProbability()};
        const int cls{(pShower > pTrack) ? 1 : 0};
        switch (view)
        {
            case TPC_VIEW_U:
                ++m_confusionU[0][cls];
                break;
            case TPC_VIEW_V:
                ++m_confusionV[0][cls];
                break;
            case TPC_VIEW_W:
                ++m_confusionW[0][cls];
                if (cls == 0)
                    correctHitListW.emplace_back(pCaloHit);
                else
                    incorrectHitListW.emplace_back(pCaloHit);
                break;
            default:
                break;
        }
    }

    for (const CaloHit *pCaloHit : showerHitList)
    {
        HitType view{pCaloHit->GetHitType()};
        const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        const float pTrack{pLArCaloHit->GetTrackProbability()};
        const float pShower{pLArCaloHit->GetShowerProbability()};
        const int cls{(pShower > pTrack) ? 1 : 0};
        switch (view)
        {
            case TPC_VIEW_U:
                ++m_confusionU[1][cls];
                break;
            case TPC_VIEW_V:
                ++m_confusionV[1][cls];
                break;
            case TPC_VIEW_W:
                ++m_confusionW[1][cls];
                if (cls == 0)
                    incorrectHitListW.emplace_back(pCaloHit);
                else
                    correctHitListW.emplace_back(pCaloHit);
                break;
            default:
                break;
        }
    }

    if (m_visualize)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &correctHitListW, "Correct", BLACK));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &incorrectHitListW, "Incorrect", RED));

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitTruthTaggingAlgorithm::VisualizeTruthTagging(const CaloHitList &trackHitList, const CaloHitList &showerHitList) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    CaloHitList trackHitListW, showerHitListW;
    for (const CaloHit *pCaloHit : trackHitList)
        if (pCaloHit->GetHitType() == TPC_VIEW_W)
            trackHitListW.emplace_back(pCaloHit);

    for (const CaloHit *pCaloHit : showerHitList)
        if (pCaloHit->GetHitType() == TPC_VIEW_W)
            showerHitListW.emplace_back(pCaloHit);
    
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHitListW, "Track Hits", BLUE));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHitListW, "Shower Hits", RED));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitTruthTaggingAlgorithm::Run()
{
    CaloHitList trackHitList, showerHitList;

    this->TagNonMuonClusters(trackHitList, showerHitList);
    this->TagMuonClusters(trackHitList, showerHitList);

    if (m_visualize)
        this->VisualizeTruthTagging(trackHitList, showerHitList);
    if (m_writeFeatures)
        this->OutputTruthTagging(trackHitList, showerHitList);

    if (m_validate)
        this->ValidateInference(trackHitList, showerHitList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitTruthTaggingAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "MuonClusterListNames", m_muonClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "NonMuonClusterListNames", m_nonMuonClusterListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFeatures", m_writeFeatures));
    if (m_writeFeatures)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFileName", m_outputFileName));
    }


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
