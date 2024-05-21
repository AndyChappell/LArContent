/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/VertexHitAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/VertexHitAssociationAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

VertexHitAssociationAlgorithm::VertexHitAssociationAlgorithm() :
    m_minClusterLayers(4),
    m_maxGapLayers(7),
    m_fitLayers(30),
    m_maxGapDistanceSquared(10.f),
    m_minCosRelativeAngle(0.985f),
    m_maxTransverseDisplacement(2.f),
    m_maxLongitudinalDisplacement(2.f),
    m_hitSizeZ(0.3f),
    m_hitSizeX(0.5f),
    m_view(TPC_3D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexHitAssociationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    ClusterList selectedClusters;
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetOrderedCaloHitList().size() > 5)
            selectedClusters.emplace_back(pCluster);
    }

    // Get the hits and the clusters
    // Project the cluster end direction to see how close it gets to various hits
    // Keep plausible candiate matches
    // Consider groups of hits (in close proximity and "matched" to some fit trajectory
    // Decide which association is best

    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));

    typedef std::tuple<const Cluster *, CaloHitVector, CaloHitVector, FloatVector, FloatVector> ClusterAssociation;
    std::vector<ClusterAssociation> clusterAssociations;
    CaloHitSet usedHits;
    for (const Cluster *const pCluster : selectedClusters)
    {
        CaloHitVector candidateInnerHits, candidateOuterHits;
        FloatVector dlInner, dtInner, dlOuter, dtOuter;
        HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        ClusterList clusterList;
        clusterList.emplace_back(pCluster);
        const LArPointingCluster pointingCluster(pCluster, 10, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
        const CartesianVector &innerDir{pointingCluster.GetInnerVertex().GetDirection() * -1.f};
        const CartesianVector &outerDir{pointingCluster.GetOuterVertex().GetDirection() * -1.f};
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            if (usedHits.find(pCaloHit) != usedHits.end())
                continue;
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            const CartesianVector &vecInner{pos - pointingCluster.GetInnerVertex().GetPosition()};
            if (vecInner.GetMagnitude() <= 10.f)
            {
                const float innerDot{innerDir.GetDotProduct(vecInner)};
                if (innerDot >= -0.5f && innerDot <= 10.f)
                {
                    const float innerCross{innerDir.GetCrossProduct(vecInner).GetMagnitude()};
                    if (innerCross < 1.f)
                    {
                        candidateInnerHits.emplace_back(pCaloHit);
                        dlInner.emplace_back(innerDot);
                        dtInner.emplace_back(innerCross);
                        usedHits.insert(pCaloHit);
                    }
                }
            }

            const CartesianVector &vecOuter{pos - pointingCluster.GetOuterVertex().GetPosition()};
            if (vecOuter.GetMagnitude() <= 10.f)
            {
                const float outerDot{outerDir.GetDotProduct(vecOuter)};
                if (outerDot >= -0.5f && outerDot <= 10.f)
                {
                    const float outerCross{outerDir.GetCrossProduct(vecOuter).GetMagnitude()};
                    if (outerCross < 1.f)
                    {
                        candidateOuterHits.emplace_back(pCaloHit);
                        dlOuter.emplace_back(outerDot);
                        dtOuter.emplace_back(outerCross);
                        usedHits.insert(pCaloHit);
                    }
                }
            }
        }

        // sort the cluster hits along the respective inner and outer pointing cluster directions
        CaloHitList tempHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(tempHitList);
        CaloHitVector clusterHitsInner(tempHitList.begin(), tempHitList.end());
        FloatVector clusterDots;
        for (const CaloHit *pCaloHit : clusterHitsInner)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            const CartesianVector &vec{pos - pointingCluster.GetInnerVertex().GetPosition()};
            clusterDots.emplace_back(innerDir.GetDotProduct(vec));
        }
        std::vector<size_t> sortIdx(clusterDots.size());
        std::iota(sortIdx.begin(), sortIdx.end(), 0);
        std::stable_sort(sortIdx.begin(), sortIdx.end(), [&clusterDots](size_t i, size_t j) { return clusterDots[i] < clusterDots[j]; });
        this->SortByIndex<const CaloHit*>(sortIdx, clusterHitsInner);

        std::vector<size_t> sortIdx1(dlInner.size());
        std::iota(sortIdx1.begin(), sortIdx1.end(), 0);
        std::stable_sort(sortIdx1.begin(), sortIdx1.end(), [&dlInner, &dtInner](size_t i, size_t j)
            {
                if (dlInner[i] != dlInner[j])
                    return dlInner[i] < dlInner[j];
                else
                    return dtInner[i] < dtInner[j];
            });
        this->SortByIndex<const CaloHit*>(sortIdx1, candidateInnerHits);
        this->SortByIndex<float>(sortIdx1, dlInner);
        this->SortByIndex<float>(sortIdx1, dtInner);

        std::vector<size_t> sortIdx2(dlOuter.size());
        std::iota(sortIdx2.begin(), sortIdx2.end(), 0);
        std::stable_sort(sortIdx2.begin(), sortIdx2.end(), [&dlOuter, &dtOuter](size_t i, size_t j)
            {
                if (dlOuter[i] != dlOuter[j])
                    return dlOuter[i] < dlOuter[j];
                else
                    return dtOuter[i] < dtOuter[j];
            });

        this->SortByIndex<const CaloHit*>(sortIdx2, candidateOuterHits);
        this->SortByIndex<float>(sortIdx2, dlOuter);
        this->SortByIndex<float>(sortIdx2, dtOuter);

        CaloHitVector clusterHitsOuter(clusterHitsInner.rbegin(), clusterHitsInner.rend());

        if (!candidateInnerHits.empty())
        {
            clusterAssociations.emplace_back(std::make_tuple(pCluster, clusterHitsInner, candidateInnerHits, dlInner, dtInner));
        }
        if (!candidateOuterHits.empty())
        {
            clusterAssociations.emplace_back(std::make_tuple(pCluster, clusterHitsOuter, candidateOuterHits, dlOuter, dtOuter));
        }
    }

    std::vector<CaloHitList> clusterHitsList;
    for (const ClusterAssociation &assoc : clusterAssociations)
    {
        const CaloHitVector &hitVector{std::get<2>(assoc)};
        const FloatVector &dlVector{std::get<3>(assoc)};
        const FloatVector &dtVector{std::get<4>(assoc)};
        CaloHitList candidateCluster;
        int ip{-1};
        for (size_t i = 0; i < hitVector.size(); ++i)
        {
            if (ip > -1)
            {
                // Disallow merge parallel hits, avoid jumping too far, and keep tight to the principal axis
                if (((dlVector[i] - dlVector[ip]) > std::numeric_limits<float>::epsilon()) && ((dlVector[i] - dlVector[ip]) < 2.f) &&
                    (std::abs(dtVector[i] - dtVector[ip]) < 0.5f))
                {
                    const CaloHit *const pPrevCaloHit{hitVector[ip]};
                    const CaloHit *const pThisCaloHit{hitVector[i]};
                    if ((pThisCaloHit->GetPositionVector() - pPrevCaloHit->GetPositionVector()).GetMagnitude() <= 2.f)
                    {
                        candidateCluster.emplace_back(pThisCaloHit);
                        ip = i;
                    }
                }
            }
            else
            {
                bool goodStart{false};
                if (dlVector[i] < 0.f && dtVector[i] < 0.2f)
                    goodStart = true;
                else if (dlVector[i] > 0.f && dtVector[i] < 0.7f)
                    goodStart = true;
                if (goodStart)
                {
                    const CaloHit *const pThisCaloHit{hitVector[i]};
                    candidateCluster.emplace_back(pThisCaloHit);
                    ip = i;
                }
            }
        }
        if (!candidateCluster.empty())
        {
            clusterHitsList.emplace_back(candidateCluster);
        }

        const CaloHitVector clusterHits{std::get<1>(assoc)};
        const CaloHitList tempClusterHits(clusterHits.begin(), clusterHits.end());
        const CaloHitList tempNewHits(candidateCluster.begin(), candidateCluster.end());
        /*if (clusterHits.front()->GetHitType() == TPC_VIEW_W)
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &tempClusterHits, "base", BLACK));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &tempNewHits, "new", AUTOITER));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }*/
    }

    const ClusterList *pOutputClusterList{nullptr};
    std::string originalClusterListName, tempListName{"temp"};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputClusterList, tempListName));
    for (CaloHitList &caloHitList : clusterHitsList)
    {
        const Cluster *pCluster{nullptr};
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList.emplace_back(caloHitList.front());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

        caloHitList.pop_front();
        for (const CaloHit *const pAssociatedCaloHit : caloHitList)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pAssociatedCaloHit));
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, originalClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalClusterListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void VertexHitAssociationAlgorithm::SortByIndex(const std::vector<size_t> &indices, std::vector<T> &vec) const
{
    std::vector<std::pair<T, int>> temp(indices.size());
    for (size_t i = 0; i < temp.size(); ++i)
    {
        temp[i].first = vec[i];
        temp[i].second = indices[i];
    }

    std::sort(temp.begin(), temp.end(), [](const std::pair<T, int> &i, const std::pair<T, int> &j){ return i.second < j.second; });
    for (size_t i = 0; i < temp.size(); ++i)
    {
        vec[i] = temp[i].first;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexHitAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapLayers", m_maxGapLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FitLayers", m_fitLayers));

    float maxGapDistance = std::sqrt(m_maxGapDistanceSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapDistance", maxGapDistance));
    m_maxGapDistanceSquared = maxGapDistance * maxGapDistance;

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitSizeZ", m_hitSizeZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitSizeX", m_hitSizeX));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
