/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/VertexClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/VertexClusterCreationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include <Eigen/Dense>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

VertexClusterCreationAlgorithm::VertexClusterCreationAlgorithm() :
    m_caloHitListName{""},
    m_vertexListName{""},
    m_hitRadii(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexClusterCreationAlgorithm::Run()
{
    // We can't guarantee that secondary vertices exist, so allow for unitialized lists
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const VertexList *pVertexList{nullptr};
    StatusCode status{PandoraContentApi::GetList(*this, m_vertexListName, pVertexList)};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, status);
    if (status == STATUS_CODE_NOT_INITIALIZED)
        return STATUS_CODE_SUCCESS;

    if (!pCaloHitList || pCaloHitList->empty() || !pVertexList || pVertexList->empty())
        return STATUS_CODE_SUCCESS;

    for (const Vertex *const pVertex : *pVertexList)
        this->IdentifyConnectingPaths(*pCaloHitList, *pVertex);
        //this->IdentifyTrackStubs(*pCaloHitList, *pVertex);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexClusterCreationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        // Pointing cluster creation requires more than 3 hits
        if (pCluster->GetOrderedCaloHitList().size() > 3)
            clusterVector.emplace_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexClusterCreationAlgorithm::IdentifyTrackStubs(const CaloHitList &caloHitList, const Vertex &vertex) const
{
    const CartesianVector pos{LArGeometryHelper::ProjectPosition(this->GetPandora(), vertex.GetPosition(), caloHitList.front()->GetHitType())};

    // Vectorize the hits
    const CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    Eigen::MatrixXf hitMatrix(caloHitVector.size(), 2), loMatrix(caloHitVector.size(), 2), hiMatrix(caloHitVector.size(), 2);
    LArEigenHelper::Vectorize(caloHitVector, hitMatrix, loMatrix, hiMatrix);
    // Compute the distance between the hits and the vertex
    Eigen::RowVectorXf vtx(2);
    vtx << pos.GetX(), pos.GetZ();
    Eigen::MatrixXf deltas(hitMatrix.rowwise() - vtx);
    Eigen::MatrixXf norms(deltas.array().pow(2).rowwise().sum().sqrt());
    // Compute the angles of the hit extrema
    Eigen::RowVectorXf loPhis(loMatrix.rows());
    Eigen::RowVectorXf hiPhis(hiMatrix.rows());
    LArEigenHelper::GetAngles(loMatrix, vtx, loPhis);
    LArEigenHelper::GetAngles(hiMatrix, vtx, hiPhis);
    // Determine angular binning (still floating point at this stage)
    const int nBins{72};
    const float pi{static_cast<float>(M_PI)};
    loPhis = loPhis.array() * nBins / (2 * pi);
    hiPhis = hiPhis.array() * nBins / (2 * pi);

    // Allocate to bins (hits can be added to multiple bins if they're wide enough)
    std::map<int, IntVector> selectedHitBinMap;
    for (int i = 0; i < norms.rows(); ++i)
    {
        if (norms(i) <= m_hitRadii)
        {
            const int bin0{loPhis(i) < hiPhis(i) ? static_cast<int>(loPhis(i)) : static_cast<int>(hiPhis(i))};
            const int bin1{loPhis(i) < hiPhis(i) ? static_cast<int>(hiPhis(i)) : static_cast<int>(loPhis(i))};
            for (int b = bin0; b <= bin1; ++b)
                selectedHitBinMap[b].emplace_back(i);
        }
    }

    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    //PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "vtx", MAGENTA, 1));

    // Compute principal components for the candidate clusters and consider fit quality
    typedef std::tuple<float, float, float> FOM;
    typedef std::pair<int, FOM> BinFOM;
    std::vector<BinFOM> metrics;
    for (const auto &[bin, hitIds] : selectedHitBinMap)
    {
        if (hitIds.size() <= 2)
            continue;
        LArPcaHelper::WeightedPointVector pointVector;
        for (const int hitId : hitIds)
        {
            const CaloHit *const pCaloHit{caloHitVector.at(hitId)};
            pointVector.emplace_back(LArPcaHelper::WeightedPoint(pCaloHit->GetPositionVector(), 1.f / pCaloHit->GetCellSize1()));
        }
        CartesianVector centroid(0, 0, 0);
        LArPcaHelper::EigenValues eigenValues(0, 0, 0);
        LArPcaHelper::EigenVectors eigenVectors;
        LArPcaHelper::RunPca(pointVector, centroid, eigenValues, eigenVectors);
        float min{std::numeric_limits<float>::max()}, max{std::numeric_limits<float>::lowest()};
        float deviation{0.f};
        // Consider sorting the longitudinal projections and checking for gaps
        for (const auto &[point, weight] : pointVector)
        {
            const CartesianVector &dir{point - centroid};
            const float longitudinal{dir.GetDotProduct(eigenVectors[0])};
            if (longitudinal > max)
                max = longitudinal;
            if (longitudinal < min)
                min = longitudinal;
            const float transverse{dir.GetDotProduct(eigenVectors[1])};
            deviation += transverse * transverse;
            //std::cout << "dt: " << transverse << " dw: " << (1.f / weight) << std::endl;
        }
        const int nHits{static_cast<int>(pointVector.size())};
        const float length{(max - min) > 0 ? max - min : 1.f};
        const float deviationPerHit{deviation / nHits};
        const float deviationPerLength{deviation / length};
        const FOM fom{std::make_tuple(length, deviationPerHit, deviationPerLength)};
        metrics.emplace_back(std::make_pair(bin, fom));
        //std::cout << "Length: " << length << " Avg Dev: " << deviationPerHit << " " << deviationPerLength << std::endl;

        CaloHitList hits;
        for (const int hitId : hitIds)
            hits.emplace_back(caloHitVector.at(hitId));
        const CartesianVector start{centroid - eigenVectors[0] * 0.5f * length};
        const CartesianVector finish{centroid + eigenVectors[0] * 0.5f * length};
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, "bin", AUTOITER));
        //PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &finish, "pca", BLUE, 1, 1));
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    std::sort(metrics.begin(), metrics.end(), [&](const BinFOM a, const BinFOM b){ return std::get<0>(a.second) > std::get<0>(b.second); });
    IntVector duplicates;
    for (size_t i = 1; i < metrics.size(); ++i)
    {
        if (std::get<0>(metrics[i].second) == std::get<0>(metrics[i - 1].second))
        {
            IntVector hitIds1{selectedHitBinMap.at(metrics[i - 1].first)};
            IntVector hitIds2{selectedHitBinMap.at(metrics[i].first)};
            if (hitIds1.size() == hitIds2.size())
            {
                std::sort(hitIds1.begin(), hitIds1.end());
                std::sort(hitIds2.begin(), hitIds2.end());
                bool duplicate{true};
                for (size_t j = 0; j < hitIds1.size(); ++j)
                {
                    if (hitIds1[j] != hitIds2[j])
                    {
                        duplicate = false;
                        break;
                    }
                }
                if (duplicate)
                    duplicates.emplace_back(metrics[i].first);
            }
        }
    }

    for (const int id : duplicates)
    {
        for (auto iter = metrics.begin(); iter != metrics.end(); )
        {
            if ((*iter).first == id)
                iter = metrics.erase(iter);
            else
                ++iter;
        }
    }

    IntVector selectedBins, usedHitIds;
    std::map<int, IntVector> clusterHitIds;
    for (const auto &[bin, fom] : metrics)
    {
        if (std::get<1>(fom) < 0.012)
        {
            selectedBins.emplace_back(bin);
            for (const int hitId : selectedHitBinMap.at(bin))
            {
                // For now, just add the hits in turn and veto used hits
                // In future do some checking to see which cluster is the best fit for a duplicated hit
                if (std::find(usedHitIds.begin(), usedHitIds.end(), hitId) == usedHitIds.end())
                {
                    clusterHitIds[bin].emplace_back(hitId);
                    usedHitIds.emplace_back(hitId);
                }
            }
        }
    }

    /*PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "vtx", MAGENTA, 1));
    for (const auto &[bin, hitIds] : selectedHitBinMap)
    {
        CaloHitList hits;
        for (const int hitId : hitIds)
            hits.emplace_back(caloHitVector.at(hitId));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, "bin", AUTOITER));
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/

    for (const auto &[bin, hitIds] : clusterHitIds)
    {
        CaloHitVector caloHits;
        for (const int hitId : hitIds)
            caloHits.emplace_back(caloHitVector.at(hitId));
        if (caloHits.size() > 2)
            this->ClusterHits(caloHits);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexClusterCreationAlgorithm::IdentifyConnectingPaths(const CaloHitList &caloHitList, const Vertex &vertex) const
{
    const CartesianVector refAxis(1.f, 0.f, 0.f);
    const CartesianVector pos{LArGeometryHelper::ProjectPosition(this->GetPandora(), vertex.GetPosition(), caloHitList.front()->GetHitType())};
    const float tolerance{(1.5f * m_hitRadii) * (1.5f * m_hitRadii)};

    const ClusterList *pClusterList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    ClusterVector cleanClusters, innerClusters, outerClusters;
    this->GetListOfCleanClusters(pClusterList, cleanClusters);
    for (const Cluster *const pCluster : cleanClusters)
    {
        try
        {
            const LArPointingCluster pointingCluster(pCluster, 10, LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
            const CartesianVector &innerPos{pointingCluster.GetInnerVertex().GetPosition()};
            const CartesianVector &outerPos{pointingCluster.GetOuterVertex().GetPosition()};
            const float innerDistance{innerPos.GetDistanceSquared(pos)};
            const float outerDistance{outerPos.GetDistanceSquared(pos)};
            if ((innerDistance < outerDistance) && (innerPos.GetDistanceSquared(pos) <= tolerance))
                innerClusters.emplace_back(pCluster);
            else if ((outerDistance <= innerDistance) && (outerPos.GetDistanceSquared(pos) <= tolerance))
                outerClusters.emplace_back(pCluster);
        }
        catch (StatusCodeException &)
        {
        }
    }

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    for (const Cluster *const pCluster : innerClusters)
    {
        const float wirePitch{LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster))};
        const LArPointingCluster initPointingCluster(pCluster, 10, wirePitch);
        // Pointing clusters always point towards their own hits
        const CartesianVector &dir{initPointingCluster.GetInnerVertex().GetDirection() * -1.f};
        const float dot{dir.GetDotProduct(refAxis)};
        // If cluster appears relatively transverse, use larger fit window to accommodate wide hits
        const LArPointingCluster pointingCluster{dot < 0.707f ? initPointingCluster : LArPointingCluster(pCluster, 20, wirePitch)};
        const CartesianVector &clusterPos{pointingCluster.GetInnerVertex().GetPosition()};
        // Pointing clusters always point towards their own hits
        const CartesianVector &clusterDir{pointingCluster.GetInnerVertex().GetDirection() * -1.f};
        const CartesianVector &unitClusterDir{clusterDir.GetUnitVector()};
        const CartesianVector &tempVertexDir{pos - clusterPos};
        // Extend the direction vector just beyond the vertex hit centre
        const CartesianVector &vertexDir{tempVertexDir * ((tempVertexDir.GetMagnitude() +  0.5f * wirePitch) / (tempVertexDir.GetMagnitude())) };
        const CartesianVector &unitVertexDir{vertexDir.GetUnitVector()};
        const CartesianVector unitPerpendicular{-unitVertexDir.GetZ(), 0, unitVertexDir.GetX()};

        // Check that the vertex dir and cluster dir are somewhat consistent
        const float clusterVertexDot{unitClusterDir.GetDotProduct(unitVertexDir)};
        if (clusterVertexDot < 0.966f)
            continue;
        CaloHitList intersectedHits;
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector hitOffset{pCaloHit->GetPositionVector() - clusterPos};
            const float hitDot{hitOffset.GetDotProduct(unitVertexDir)};
            if ((hitDot >= 0.f) && (hitDot < vertexDir.GetMagnitude()))
            {
                const CartesianVector llOffset{(pCaloHit->GetPositionVector() + CartesianVector(-pCaloHit->GetCellSize1() * 0.5f, 0.f, -0.5f * wirePitch)) - clusterPos};
                const CartesianVector hhOffset{(pCaloHit->GetPositionVector() + CartesianVector(+pCaloHit->GetCellSize1() * 0.5f, 0.f, +0.5f * wirePitch)) - clusterPos};
                const float llDot{llOffset.GetDotProduct(unitPerpendicular)};
                const float hhDot{hhOffset.GetDotProduct(unitPerpendicular)};
                if (llDot * hhDot <= 0.f)
                {
                    intersectedHits.emplace_back(pCaloHit);
                }
                else
                {
                    const CartesianVector lhOffset{(pCaloHit->GetPositionVector() + CartesianVector(-pCaloHit->GetCellSize1() * 0.5f, 0.f, +0.5f * wirePitch)) - clusterPos};
                    const CartesianVector hlOffset{(pCaloHit->GetPositionVector() + CartesianVector(+pCaloHit->GetCellSize1() * 0.5f, 0.f, -0.5f * wirePitch)) - clusterPos};
                    const float lhDot{lhOffset.GetDotProduct(unitPerpendicular)};
                    const float hlDot{hlOffset.GetDotProduct(unitPerpendicular)};
                    if (lhDot * hlDot <= 0.f)
                    {
                        intersectedHits.emplace_back(pCaloHit);
                    }
                }
            }
        }
        // Check for continuity of hits, both within the selected hits and with respect to the base cluster
        // The hits should not start a large distance from the base cluster, and should not contain large gaps

        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "vtx", MAGENTA, 1));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "hits", GRAY));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &intersectedHits, "path_hits", ORANGE));
        ClusterList clusters({pCluster});
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusters, "cluster", RED));
        const CartesianVector a(clusterPos + vertexDir);
        const CartesianVector b(clusterPos + clusterDir);
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &clusterPos, &a, "vertex_dir", BLACK, 1, 1));
        PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &clusterPos, &b, "cluster_dir", BLUE, 1, 1));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

/*    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "vtx", MAGENTA, 1));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "hits", GRAY));
    ClusterList nearbyInnerClusters(innerClusters.begin(), innerClusters.end());
    ClusterList nearbyOuterClusters(outerClusters.begin(), outerClusters.end());
    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &nearbyInnerClusters, "inner", BLUE));
    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &nearbyOuterClusters, "outer", RED));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexClusterCreationAlgorithm::ClusterHits(const CaloHitVector &caloHitVector) const
{
    const Cluster *pCluster{nullptr};
    const ClusterList *pOutputClusterList{nullptr};
    std::string originalClusterListName, tempListName{"temp"};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalClusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputClusterList, tempListName));

    for (const CaloHit *const pCaloHit : caloHitVector)
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        if (!pCluster)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.emplace_back(pCaloHit);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit));
        }
    }
    if (pCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, originalClusterListName));
    }
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalClusterListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexClusterCreationAlgorithm::GetClosestApproach(const CaloHitVector &caloHitVector, const CartesianVector &vertex) const
{
    float closestApproach{std::numeric_limits<float>::max()};
    for (const CaloHit *const pCaloHit : caloHitVector)
    {
        const CartesianVector delta{pCaloHit->GetPositionVector() - vertex};
        float dr{delta.GetMagnitudeSquared()};
        if (dr < closestApproach)
            closestApproach = dr;
    }

    return closestApproach;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitRadii", m_hitRadii));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
