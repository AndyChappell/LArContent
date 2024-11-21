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
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const VertexList *pVertexList{nullptr};
    StatusCode status{PandoraContentApi::GetList(*this, m_vertexListName, pVertexList)};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, status);
    if (status == STATUS_CODE_NOT_INITIALIZED)
        return STATUS_CODE_SUCCESS;

    if (!pCaloHitList || pCaloHitList->empty() || !pVertexList || pVertexList->empty())
        return STATUS_CODE_SUCCESS;

    for (const Vertex *const pVertex : *pVertexList)
        this->IdentifyTrackStubs(*pCaloHitList, *pVertex);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexClusterCreationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetOrderedCaloHitList().size() > 2)
            clusterVector.emplace_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexClusterCreationAlgorithm::IdentifyTrackStubs(const CaloHitList &caloHitList, const Vertex &vertex) const
{
    const CartesianVector pos{LArGeometryHelper::ProjectPosition(this->GetPandora(), vertex.GetPosition(), caloHitList.front()->GetHitType())};

    const ClusterList *pClusterList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    ClusterVector cleanClusters, innerClusters, outerClusters;
    this->GetListOfCleanClusters(pClusterList, cleanClusters);
    for (const Cluster *const pCluster : cleanClusters)
    {
        const LArPointingCluster pointingCluster(pCluster, 10, LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
        const CartesianVector &innerPos{pointingCluster.GetInnerVertex().GetPosition()};
        if (innerPos.GetDistanceSquared(pos) <= (1.5f * m_hitRadii) * (1.5f * m_hitRadii))
            innerClusters.emplace_back(pCluster);
        const CartesianVector &outerPos{pointingCluster.GetOuterVertex().GetPosition()};
        if (outerPos.GetDistanceSquared(pos) <= (1.5f * m_hitRadii) * (1.5f * m_hitRadii))
            outerClusters.emplace_back(pCluster);
    }

    // Vectorize thie hits
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
    std::map<int, std::vector<int>> selectedHitBinMap;
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

    /*bool clusterMade{true};
    while (clusterMade)
    {
        clusterMade = false;
        for (const int bin : sortedBins)
        {
            int target{bin};
            int binAlt1{bin > 0 ? bin - 1 : nBins - 1};
            int binAlt2{bin < (nBins - 1) ? bin + 1 : 0};
            float delta{closestApproach[bin]};
            float deltaAlt{closestApproach[binAlt1]};
            if (deltaAlt < delta)
            {
                target = binAlt1;
                delta = deltaAlt;
            }
            deltaAlt = closestApproach[binAlt2];
            if (deltaAlt < delta)
            {
                target = binAlt2;
                delta = deltaAlt;
            }
            // Candidate clusters are sorted by the number of hits, if there aren't two hits there's nothing to cluster
            if (hitBins[target] < 2)
            {
                hitBins[target] = 0;
                closestApproach[target] = std::numeric_limits<float>::max();
                // Need to allow for sparse adjacent bins that happen to be closer than the nominal bin
                if (target != bin)
                    clusterMade = true;
                std::sort(sortedBins.begin(), sortedBins.end(), SortBinsFunc);
                break;
            }
            const CaloHit *const pSeedHit{caloHits[target].front()};
            caloHits[target].erase(caloHits[target].begin());
            clusterMade = this->ClusterHits(pSeedHit, caloHits[target], pos);
            if (clusterMade)
            {
                caloHits[target].erase(std::remove_if(caloHits[target].begin(), caloHits[target].end(),
                    [this](const CaloHit *const pCaloHit)
                    {
                        return !PandoraContentApi::IsAvailable(*this, pCaloHit);
                    }), caloHits[target].end());
                hitBins[target] = caloHits[target].size();
                closestApproach[target] = this->GetClosestApproach(caloHits[target], pos);
                std::sort(sortedBins.begin(), sortedBins.end(), SortBinsFunc);
                break;
            }
            else
            {
                // If we didn't make a cluster re-insert the seed into the hit list
                caloHits[target].insert(caloHits[target].begin(), pSeedHit);
            }
        }
    }*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexClusterCreationAlgorithm::ClusterHits(const CaloHit *const pSeedHit, const CaloHitVector &caloHitVector,
    const CartesianVector &vertex) const
{
    // We know the seed hit is not null and the hit vector must have at least one element
    CaloHitVector clusteredHits;
    clusteredHits.emplace_back(pSeedHit);
    CartesianVector clusterDirection{(pSeedHit->GetPositionVector() - vertex).GetUnitVector()};

    for (size_t i = 0; i < caloHitVector.size(); ++i)
    {
        const CaloHit *pCurrentHit{caloHitVector.at(i)};
        const int N{static_cast<int>(clusteredHits.size())};
        const CaloHit *pAnchorHit{nullptr};
        float tolerance{0.9969f};
        switch (N)
        {
            case 1:
                pAnchorHit = clusteredHits[0];
                tolerance = 0.985f;
                break;
            case 2:
                pAnchorHit = clusteredHits[0];
                tolerance = 0.9875f;
                break;
            case 3:
                pAnchorHit = clusteredHits[0];
                tolerance = 0.99f;
                break;
            case 4:
                pAnchorHit = clusteredHits[0];
                tolerance = 0.9925f;
                break;
            default:
                pAnchorHit = clusteredHits[N - 5];
                tolerance = 0.9925f;
                break;
        }
        CartesianVector localDirection{(pCurrentHit->GetPositionVector() - pAnchorHit->GetPositionVector()).GetUnitVector()};
        const float dr{(pCurrentHit->GetPositionVector() - clusteredHits[N - 1]->GetPositionVector()).GetMagnitudeSquared()};

        if (dr <= 4.f)
        {
            // If we aren't too far from the last cluster and the direction is within tight angular constraints, cluster
            // Also allow some latitude if we've only clustered the seed hit so far
            const float dotProduct{clusterDirection.GetDotProduct(localDirection)};
            if (dotProduct >= tolerance)
                clusteredHits.emplace_back(pCurrentHit);
        }
        else
        {
            // No downstream hits can be closer than the last, we're done
            break;
        }
    }

    const Cluster *pCluster{nullptr};
    if (clusteredHits.size() > 1)
    {
        const ClusterList *pOutputClusterList{nullptr};
        std::string originalClusterListName, tempListName{"temp"};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalClusterListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputClusterList, tempListName));

        for (const CaloHit *const pCaloHit : clusteredHits)
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

    return pCluster;
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
