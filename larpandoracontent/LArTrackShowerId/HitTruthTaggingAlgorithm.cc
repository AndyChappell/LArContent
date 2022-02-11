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
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTrackShowerId/HitTruthTaggingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

HitTruthTaggingAlgorithm::HitTruthTaggingAlgorithm() : m_minTrackRatio{0.5f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitTruthTaggingAlgorithm::TagMuonClusters() const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    CaloHitList trackHitList, showerHitList, diffuseHitList;

    const int SHOWER{1}, TRACK{2}, DIFFUSE{4};
    for (std::string clusterListName : m_muonClusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

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
                try
                {
                    const CaloHit *pCaloHit{caloHitList.front()};
                    const MCParticle *const pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
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
                catch (const StatusCodeException &)
                {
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
                try
                {
                    const CaloHit *pCaloHit{caloHitList.front()};
                    const MCParticle *const pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                    const LArMCParticle *const pLArMCParticle{dynamic_cast<const LArMCParticle *>(pMCParticle)};
                    if (!pLArMCParticle)
                        continue;

                    if (LArMCParticleHelper::IsCapture(pLArMCParticle) || LArMCParticleHelper::IsNuclear(pLArMCParticle) ||
                        LArMCParticleHelper::IsIonisation(pLArMCParticle))
                    {
                        tag = DIFFUSE;
                    }
                    else if(LArMCParticleHelper::IsInelasticScatter(pLArMCParticle))
                    {
                        if (caloHitList.size() >= std::max(5.f, m_minTrackRatio * mcToLargestClusterMap[pMCParticle]))
                            tag = TRACK;
                        else
                            tag = DIFFUSE;
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
                    if (pCaloHit->GetHitType() == TPC_VIEW_W)
                    {
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
                                case DIFFUSE:
                                    diffuseHitList.emplace_back(pHit);
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                }
                catch (const StatusCodeException &)
                {
                }
            }
        }
    }
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHitList, "Track Hits", BLUE));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHitList, "Shower Hits", RED));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &diffuseHitList, "Diffuse Hits", MAGENTA));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitTruthTaggingAlgorithm::TagNonMuonClusters() const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    CaloHitList trackHitList, showerHitList, michelHitList, diffuseHitList;

    const int SHOWER{1}, TRACK{2}, MICHEL{3}, DIFFUSE{4};
    for (std::string clusterListName : m_nonMuonClusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

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
                        tag = DIFFUSE;
                    }
                    else if(LArMCParticleHelper::IsInelasticScatter(pLArMCParticle))
                    {
                        if (caloHitList.size() >= 5)
                        {
                            if (pdg == PHOTON || pdg == E_MINUS)
                                tag = DIFFUSE;
                            else
                                tag = TRACK;
                        }
                        else
                        {
                            tag = DIFFUSE;
                        }
                    }
                    else if (pdg == PHOTON)
                    {
                        tag = SHOWER;
                    }
                    else if (pdg == E_MINUS)
                    {
                        if (LArMCParticleHelper::IsDecay(pLArMCParticle) && std::abs(pLArMCParticle->GetParentList().front()->GetParticleId()) == MU_MINUS)
                            tag = MICHEL;
                        else
                            tag = SHOWER;
                    }
                    else
                    {
                        tag = TRACK;
                    }
                    if (pCaloHit->GetHitType() == TPC_VIEW_W)
                    {
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
                                case MICHEL:
                                    michelHitList.emplace_back(pHit);
                                    break;
                                case DIFFUSE:
                                    diffuseHitList.emplace_back(pHit);
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                }
                catch (const StatusCodeException &)
                {
                }
            }
        }
    }

    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &trackHitList, "Track Hits", BLUE));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHitList, "Shower Hits", RED));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &michelHitList, "Michel Hits", BLUE));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &diffuseHitList, "Diffuse Hits", MAGENTA));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitTruthTaggingAlgorithm::Run()
{
    this->TagNonMuonClusters();
    this->TagMuonClusters();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitTruthTaggingAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "MuonClusterListNames", m_muonClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "NonMuonClusterListNames", m_nonMuonClusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
