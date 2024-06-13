/**
 *  @file   larpandoracontent/LArTwoDReco/LArAssociatedCluster/DlVertexAssociatedClusterAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoradlcontent/LArVertex/DlVertexAssociatedClusterAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlVertexAssociatedClusterAlgorithm::DlVertexAssociatedClusterAlgorithm() :
    m_clusterListName{""},
    m_vertexListName{""},
    m_rInfluenceSquared{25.f},
    m_halfBoundingEdge{16.f},
    m_nBins{64},
    m_trainingMode{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlVertexAssociatedClusterAlgorithm::~DlVertexAssociatedClusterAlgorithm()
{
    if (m_trainingMode)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlVertexAssociatedClusterAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexAssociatedClusterAlgorithm::Run()
{
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));
    const VertexList *pVertexList{nullptr};
    StatusCode status{PandoraContentApi::GetList(*this, m_vertexListName, pVertexList)};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, status);

    if (status == STATUS_CODE_NOT_INITIALIZED)
        return STATUS_CODE_SUCCESS;

    if (!pClusterList || pClusterList->empty() || !pVertexList || pVertexList->empty())
        return STATUS_CODE_SUCCESS;

    this->IdentifyAssociatedClusters(*pClusterList, *pVertexList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DlVertexAssociatedClusterAlgorithm::IdentifyAssociatedClusters(const ClusterList &clusterList, const VertexList &vertexList) const
{
    typedef std::map<const Cluster *, MCParticleList> ClusterMCMap;
    std::map<const Cluster *, CaloHitList> clusterToHitsMap;
    ClusterMCMap clusterToMCMap;
    for (const Cluster *pCluster : clusterList)
    {
        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);
        if (clusterHits.empty())
            continue;
        clusterToHitsMap[pCluster] = clusterHits;
        MCParticleWeightMap clusterMCContrib;
        double totalWeight{0.f};
        for (const CaloHit *const pCaloHit : clusterHits)
        {
            const MCParticleWeightMap &hitMCContrib{pCaloHit->GetMCParticleWeightMap()};
            for (const auto &[pMC, weight] : hitMCContrib)
            {
                if (clusterMCContrib.find(pMC) == clusterMCContrib.end())
                    clusterMCContrib[pMC] = 0.;
                clusterMCContrib[pMC] += weight;
                totalWeight += weight;
            }
        }
        for (const auto &[pMC, weight] : clusterMCContrib)
        {
            if (weight > 0.3 * totalWeight)
                clusterToMCMap[pCluster].emplace_back(pMC);
        }
    }
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};

    typedef std::map<const Cluster *, CaloHitList> ClusterHitMap;
    // Gather the hits within a given radii of the vertices
    for (const Vertex *const pVertex : vertexList)
    {
        ClusterHitMap vertexHits, activeHits, backgroundHits;
        const CartesianVector &pos3D{pVertex->GetPosition()};
        CartesianVector pos(0, 0, 0), low(0, 0, 0), high(0, 0, 0);
        for (const Cluster *pCluster : clusterList)
        {
            const CaloHitList &clusterHits{clusterToHitsMap.at(pCluster)};
            HitType view{clusterHits.front()->GetHitType()};
            switch (view)
            {
                case TPC_VIEW_U:
                    pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoU(pos3D.GetY(), pos3D.GetZ())));
                    break;
                case TPC_VIEW_V:
                    pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoV(pos3D.GetY(), pos3D.GetZ())));
                    break;
                case TPC_VIEW_W:
                    pos.SetValues(pos3D.GetX(), 0.f, static_cast<float>(transform->YZtoW(pos3D.GetY(), pos3D.GetZ())));
                    break;
                default:
                    return;
            }
            low.SetValues(pos.GetX() - m_halfBoundingEdge, 0.f, pos.GetZ() - m_halfBoundingEdge);
            high.SetValues(pos.GetX() + m_halfBoundingEdge, 0.f, pos.GetZ() + m_halfBoundingEdge);

            for (const CaloHit *const pCaloHit : clusterToHitsMap.at(pCluster))
            {
                const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
                const float distanceSquared{hitPos.GetDistanceSquared(pos)};
                if (distanceSquared <= m_rInfluenceSquared)
                {
                    vertexHits[pCluster].emplace_back(pCaloHit);
                }
                else
                {
                    if (low.GetX() < hitPos.GetX() && low.GetZ() < hitPos.GetZ() && hitPos.GetX() < high.GetX() && hitPos.GetZ() < high.GetZ())
                        activeHits[pCluster].emplace_back(pCaloHit);
                }
            }

            for (const Cluster *pOtherCluster : clusterList)
            {
                if (pOtherCluster == pCluster)
                    continue;
                for (const CaloHit *const pCaloHit : clusterToHitsMap.at(pOtherCluster))
                {
                    const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
                    if (low.GetX() < hitPos.GetX() && low.GetZ() < hitPos.GetZ() && hitPos.GetX() < high.GetX() && hitPos.GetZ() < high.GetZ())
                        backgroundHits[pCluster].emplace_back(pCaloHit);
                }
            }
        }
        int id{0};
        std::vector<float> imageVertex(m_nBins * m_nBins, 0.f);
        std::vector<float> imageActive(m_nBins * m_nBins, 0.f);
        std::vector<float> imageBackground(m_nBins * m_nBins, 0.f);
        std::vector<float> cls(m_nBins * m_nBins, 0.f);
        for (const Cluster *pCluster : clusterList)
        {
            if (!vertexHits[pCluster].empty() && !activeHits[pCluster].empty())
            {
                imageVertex.assign(m_nBins * m_nBins, false);
                imageActive.assign(m_nBins * m_nBins, false);
                imageBackground.assign(m_nBins * m_nBins, false);
                cls.assign(m_nBins * m_nBins, false);
                ++id;
                for (const CaloHit *const pCaloHit : vertexHits[pCluster])
                {
                    bool isTarget{false};
                    const MCParticleWeightMap &hitMCContrib{pCaloHit->GetMCParticleWeightMap()};
                    double totalWeight{0.};
                    for (const auto &[pMC, weight] : hitMCContrib)
                        totalWeight += weight;
                    for (const auto &[pMC, weight] : hitMCContrib)
                    {
                        if (weight > 0.3 * totalWeight)
                        {
                            const MCParticleList &mcList{clusterToMCMap[pCluster]};
                            if (std::find(mcList.begin(), mcList.end(), pMC) != mcList.end())
                            {
                                isTarget = true;
                                break;
                            }
                        }
                    }
                    const int row{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetZ() - low.GetZ()) / (high.GetZ() - low.GetZ()))};
                    const int col{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetX() - low.GetX()) / (high.GetX() - low.GetX()))};
                    imageVertex[row * m_nBins + col] += pCaloHit->GetMipEquivalentEnergy();
                    cls[row * m_nBins + col] += isTarget ? pCaloHit->GetMipEquivalentEnergy() : -pCaloHit->GetMipEquivalentEnergy();
                }
                for (const CaloHit *const pCaloHit : activeHits[pCluster])
                {
                    bool isTarget{false};
                    const MCParticleWeightMap &hitMCContrib{pCaloHit->GetMCParticleWeightMap()};
                    double totalWeight{0.};
                    for (const auto &[pMC, weight] : hitMCContrib)
                        totalWeight += weight;
                    for (const auto &[pMC, weight] : hitMCContrib)
                    {
                        if (weight > 0.3 * totalWeight)
                        {
                            const MCParticleList &mcList{clusterToMCMap[pCluster]};
                            if (std::find(mcList.begin(), mcList.end(), pMC) != mcList.end())
                            {
                                isTarget = true;
                                break;
                            }
                        }
                    }
                    const int row{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetZ() - low.GetZ()) / (high.GetZ() - low.GetZ()))};
                    const int col{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetX() - low.GetX()) / (high.GetX() - low.GetX()))};
                    imageActive[row * m_nBins + col] += pCaloHit->GetMipEquivalentEnergy();
                    cls[row * m_nBins + col] += isTarget ? pCaloHit->GetMipEquivalentEnergy() : -pCaloHit->GetMipEquivalentEnergy();
                }
                for (const CaloHit *const pCaloHit : backgroundHits[pCluster])
                {
                    bool isTarget{false};
                    const MCParticleWeightMap &hitMCContrib{pCaloHit->GetMCParticleWeightMap()};
                    double totalWeight{0.};
                    for (const auto &[pMC, weight] : hitMCContrib)
                        totalWeight += weight;
                    for (const auto &[pMC, weight] : hitMCContrib)
                    {
                        if (weight > 0.3 * totalWeight)
                        {
                            const MCParticleList &mcList{clusterToMCMap[pCluster]};
                            if (std::find(mcList.begin(), mcList.end(), pMC) != mcList.end())
                            {
                                isTarget = true;
                                break;
                            }
                        }
                    }
                    const int row{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetZ() - low.GetZ()) / (high.GetZ() - low.GetZ()))};
                    const int col{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetX() - low.GetX()) / (high.GetX() - low.GetX()))};
                    imageBackground[row * m_nBins + col] += pCaloHit->GetMipEquivalentEnergy();
                    cls[row * m_nBins + col] += isTarget ? pCaloHit->GetMipEquivalentEnergy() : -pCaloHit->GetMipEquivalentEnergy();
                }

                /*PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos3D, "true", BLACK, 3));
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vertexHits[pCluster], "vertex", RED));
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &activeHits[pCluster], "active", BLUE));
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &backgroundHits[pCluster], "background", GRAY));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
            }

            IntVector indices, classes;
            FloatVector vertexWeights, activeWeights, backgroundWeights;
            for (int index = 0; index < m_nBins * m_nBins; ++index)
            {
                const float hit{std::max({imageVertex[index], imageActive[index], imageBackground[index]})};
                if (hit)
                {
                    indices.emplace_back(index);
                    classes.emplace_back(cls[index] > 0 ? 1 : 0);
                    vertexWeights.emplace_back(imageVertex[index] > 0 ? imageVertex[index] : 0);
                    activeWeights.emplace_back(imageActive[index] > 0 ? imageActive[index] : 0);
                    backgroundWeights.emplace_back(imageBackground[index] > 0 ? imageBackground[index] : 0);
                }
            }
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "indices", &indices));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "wVtx", &vertexWeights));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "wAct", &activeWeights));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "wBak", &backgroundWeights));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_rootTreeName, "classes", &classes));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_rootTreeName));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexAssociatedClusterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RadiusOfInfluenceSquared", m_rInfluenceSquared));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HalfBoundingEdge", m_halfBoundingEdge));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NBins", m_nBins));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    }

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
