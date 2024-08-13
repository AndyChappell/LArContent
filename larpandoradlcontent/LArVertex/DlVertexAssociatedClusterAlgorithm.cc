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
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArUtilityHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

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
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexAssociatedClusterAlgorithm::Infer()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));
    const VertexList *pVertexList{nullptr};
    StatusCode status{PandoraContentApi::GetList(*this, m_vertexListName, pVertexList)};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, status);

    if (status == STATUS_CODE_NOT_INITIALIZED)
        return STATUS_CODE_SUCCESS;

    if (!pClusterList || pClusterList->empty() || !pVertexList || pVertexList->empty())
        return STATUS_CODE_SUCCESS;

    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};

    typedef std::map<const Cluster *, CaloHitList> ClusterHitMap;
    typedef std::map<const Cluster *, ClusterSet> ClusterClusterMap;
    typedef std::map<const Cluster *, int> ClusterVoteMap;
    typedef std::map<const Cluster *, CaloHitSet> ClusterToHitSetMap;

    ClusterVoteMap clusterBalance;
    ClusterToHitSetMap clusterToBadHitsMap;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));
    if (!pClusterList || pClusterList->empty())
        return STATUS_CODE_SUCCESS;

    std::map<const Cluster *, CaloHitList> clusterToHitsMap;
    for (const Cluster *pCluster : *pClusterList)
    {
        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);
        if (clusterHits.empty())
            continue;
        clusterToHitsMap[pCluster] = clusterHits;
    }

    // Gather the hits within a given radii of the vertices
    for (const Vertex *const pVertex : *pVertexList)
    {
        ClusterHitMap targetHits, backgroundHits;
        ClusterClusterMap contextClusters;
        const CartesianVector &pos3D{pVertex->GetPosition()};
        CartesianVector pos(0, 0, 0), low(0, 0, 0), high(0, 0, 0);
        for (const Cluster *pCluster : *pClusterList)
        {
            if (clusterToHitsMap.find(pCluster) == clusterToHitsMap.end())
                continue;
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
                    continue;
            }
            low.SetValues(pos.GetX() - m_halfBoundingEdge, 0.f, pos.GetZ() - m_halfBoundingEdge);
            high.SetValues(pos.GetX() + m_halfBoundingEdge, 0.f, pos.GetZ() + m_halfBoundingEdge);

            for (const CaloHit *const pCaloHit : clusterToHitsMap.at(pCluster))
            {
                const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
                const float distanceSquared{hitPos.GetDistanceSquared(pos)};
                if (distanceSquared <= m_rInfluenceSquared)
                {
                    targetHits[pCluster].emplace_back(pCaloHit);
                }
                else
                {
                    if (low.GetX() < hitPos.GetX() && low.GetZ() < hitPos.GetZ() && hitPos.GetX() < high.GetX() && hitPos.GetZ() < high.GetZ())
                        targetHits[pCluster].emplace_back(pCaloHit);
                }
            }

            for (const Cluster *pOtherCluster : *pClusterList)
            {
                if (clusterToHitsMap.find(pOtherCluster) == clusterToHitsMap.end())
                    continue;
                if (pOtherCluster == pCluster)
                    continue;
                for (const CaloHit *const pCaloHit : clusterToHitsMap.at(pOtherCluster))
                {
                    const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
                    if (low.GetX() < hitPos.GetX() && low.GetZ() < hitPos.GetZ() && hitPos.GetX() < high.GetX() && hitPos.GetZ() < high.GetZ())
                    {
                        backgroundHits[pCluster].emplace_back(pCaloHit);
                        if (contextClusters[pCluster].find(pOtherCluster) == contextClusters[pCluster].end())
                            contextClusters[pCluster].insert(pOtherCluster);
                    }
                }
            }
        }

        typedef std::pair<int, int> Pixel;
        typedef std::map<const CaloHit *, Pixel> HitToPixelMap;
        HitToPixelMap hitToPixelMap;
        for (const Cluster *pCluster : *pClusterList)
        {
            if (!targetHits[pCluster].empty() && targetHits[pCluster].size() >= 5)
            {
                LArDLHelper::TorchInput networkInput;
                LArDLHelper::InitialiseInput({1, 2, m_nBins, m_nBins}, networkInput);
                auto accessor = networkInput.accessor<float, 4>();

                // Need to create a pixel map for post-processing
                for (const CaloHit *const pCaloHit : targetHits[pCluster])
                {
                    const int row{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetZ() - low.GetZ()) / (high.GetZ() - low.GetZ()))};
                    const int col{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetX() - low.GetX()) / (high.GetX() - low.GetX()))};
                    accessor[0][0][row][col] += pCaloHit->GetMipEquivalentEnergy();
                    hitToPixelMap[pCaloHit] = Pixel(row, col);
                }
                for (const CaloHit *const pCaloHit : backgroundHits[pCluster])
                {
                    const int row{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetZ() - low.GetZ()) / (high.GetZ() - low.GetZ()))};
                    const int col{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetX() - low.GetX()) / (high.GetX() - low.GetX()))};
                    accessor[0][1][row][col] += pCaloHit->GetMipEquivalentEnergy();
                    hitToPixelMap[pCaloHit] = Pixel(row, col);
                }

                LArDLHelper::TorchInputVector inputs;
                inputs.push_back(networkInput);
                LArDLHelper::TorchOutput output;
                LArDLHelper::Forward(m_model, inputs, output);

                auto classes{torch::argmax(output, 1)};
                // the argmax result is a 1 x height x width tensor where each element is a class id
                auto classesAccessor{classes.accessor<long, 3>()};
                CaloHitList targetTrue, targetFalse;
                for (const CaloHit *const pCaloHit : targetHits[pCluster])
                {
                    const int row{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetZ() - low.GetZ()) / (high.GetZ() - low.GetZ()))};
                    const int col{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetX() - low.GetX()) / (high.GetX() - low.GetX()))};
                    const auto cls{classesAccessor[0][row][col]};
                    if (cls == 1)
                    {
                        targetTrue.emplace_back(pCaloHit);
                    }
                    else
                    {
                        targetFalse.emplace_back(pCaloHit);
                        clusterToBadHitsMap[pCluster].insert(pCaloHit);
                    }
                }
                //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &targetHits[pCluster], "Cluster", BLACK));
                //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &targetTrue, "Good Hits", BLUE));
                //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &targetFalse, "Bad Hits", RED));
                //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                if (clusterBalance.find(pCluster) == clusterBalance.end())
                    clusterBalance[pCluster] = 0;
                if (targetFalse.empty())
                    ++clusterBalance[pCluster];
                else
                    --clusterBalance[pCluster];
            }
        }
    }

    for (auto &[pCluster, balance] : clusterBalance)
    {
        //std::cout << "Cluster: " << pCluster << " Balance: " << balance << std::endl;
        if (balance >= 0)
            continue;
        const CaloHitList clusterHits{clusterToHitsMap.at(pCluster)};
        CartesianVector centroid(0, 0, 0);
        LArPcaHelper::EigenValues eigenValues(0, 0, 0);
        LArPcaHelper::EigenVectors eigenVectors;
        LArPcaHelper::RunPca(clusterHits, centroid, eigenValues, eigenVectors);
        FloatVector projection, isGood;
        CaloHitVector sortedHits;
        for (const CaloHit *const pCaloHit : clusterHits)
        {
            sortedHits.emplace_back(pCaloHit);
            if (std::find(clusterToBadHitsMap[pCluster].begin(), clusterToBadHitsMap[pCluster].end(), pCaloHit) == clusterToBadHitsMap[pCluster].end())
                isGood.emplace_back(1);
            else
                isGood.emplace_back(0);
            const CartesianVector posVec{pCaloHit->GetPositionVector() - centroid};
            float dot{posVec.GetDotProduct(eigenVectors[0])};
            projection.emplace_back(dot);
        }

        auto comparator{[](const float a, const float b){ return a < b; }};
        std::vector<size_t> order{LArUtilityHelper::GetSortIndices(projection, comparator)};
        LArUtilityHelper::SortByIndices(order, sortedHits);
        LArUtilityHelper::SortByIndices(order, projection);
        LArUtilityHelper::SortByIndices(order, isGood);
        size_t currentStreak{0}, longestStreak{0};
        size_t currentStart{0}, longestStart{0};
        for (size_t i = 0; i < isGood.size(); ++i)
        {
            if (isGood[i])
            {
                ++currentStreak;
            }
            else
            {
                if (currentStreak > longestStreak)
                {
                    longestStreak = currentStreak;
                    longestStart = currentStart;
                }
                currentStart = i + 1;
                currentStreak = 0;
            }
        }
        if (currentStreak > longestStreak)
        {
            longestStreak = currentStreak;
            longestStart = currentStart;
        }
        const CaloHitSet badHits{clusterToBadHitsMap[pCluster]};
        PandoraContentApi::Cluster::Parameters firstParameters, secondParameters;
        for (size_t i = 0; i < isGood.size(); ++i)
        {
            const CaloHit *const pCaloHit{sortedHits[i]};
            if (i >= longestStart && i < (longestStart + longestStreak))
            {
                firstParameters.m_caloHitList.emplace_back(pCaloHit);
            }
            else
            {
                secondParameters.m_caloHitList.emplace_back(pCaloHit);
            }
        }
        //const CaloHitList &originalHits{clusterToHitsMap.at(pCluster)};
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &originalHits, "Original Cluster", BLACK));
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &firstParameters.m_caloHitList, "Retained Hits", BLUE));
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &secondParameters.m_caloHitList, "Split Hits", RED));
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        if (firstParameters.m_caloHitList.empty() || secondParameters.m_caloHitList.empty())
            continue;

        // Begin cluster fragmentation operations
        const ClusterList clusterList(1, pCluster);
        std::string clusterListToSaveName, clusterListToDeleteName;

        if (STATUS_CODE_SUCCESS != PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
            clusterListToSaveName))
            continue;
        // Create new clusters
        const Cluster *pFirstCluster(nullptr), *pSecondCluster(nullptr);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::Cluster::Create(*this, firstParameters, pFirstCluster))
            continue;
        clusterToHitsMap[pFirstCluster] = firstParameters.m_caloHitList;
        if (STATUS_CODE_SUCCESS != PandoraContentApi::Cluster::Create(*this, secondParameters, pSecondCluster))
            continue;
        clusterToHitsMap[pSecondCluster] = secondParameters.m_caloHitList;

        // End cluster fragmentation operations
        if (STATUS_CODE_SUCCESS != PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName))
            continue;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlVertexAssociatedClusterAlgorithm::PrepareTrainingSample() const
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
        FloatVector imageVertex(m_nBins * m_nBins, 0.f);
        FloatVector imageActive(m_nBins * m_nBins, 0.f);
        FloatVector imageBackground(m_nBins * m_nBins, 0.f);
        FloatVector cls(m_nBins * m_nBins, 0.f);
        for (const Cluster *pCluster : clusterList)
        {
            if (!vertexHits[pCluster].empty() && !activeHits[pCluster].empty())
            {
                imageVertex.assign(m_nBins * m_nBins, false);
                imageActive.assign(m_nBins * m_nBins, false);
                imageBackground.assign(m_nBins * m_nBins, false);
                cls.assign(m_nBins * m_nBins, false);
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
                    const int row{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetZ() - low.GetZ()) / (high.GetZ() - low.GetZ()))};
                    const int col{static_cast<int>(m_nBins * (pCaloHit->GetPositionVector().GetX() - low.GetX()) / (high.GetX() - low.GetX()))};
                    imageBackground[row * m_nBins + col] += pCaloHit->GetMipEquivalentEnergy();
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
                    if (cls[index] > 0)
                        classes.emplace_back(1);
                    else if (cls[index] < 0)
                        classes.emplace_back(2);
                    else
                        classes.emplace_back(0);
                    vertexWeights.emplace_back(imageVertex[index]);
                    activeWeights.emplace_back(imageActive[index]);
                    backgroundWeights.emplace_back(imageBackground[index]);
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
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_model);
    }

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
