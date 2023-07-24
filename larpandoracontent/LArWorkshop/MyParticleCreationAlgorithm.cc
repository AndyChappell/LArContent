/**
 *  @file   larpandoracontent/LArWorkshop/MyParticleCreationAlgorithm.cc
 *
 *  @brief  Implementation of a custom particle creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArWorkshop/MyParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MyParticleCreationAlgorithm::MyParticleCreationAlgorithm() : m_fitWindow{20.f}, m_nSamplingPoints{10}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyParticleCreationAlgorithm::Run()
{
    ClusterMap clusterMap;

    for (std::string clusterListName : m_clusterListNames)
    {
        ClusterVector clusterVector;
        // ATTN: Need to check the clusters have enough hits for the sliding fit
        this->GetLongClusters(clusterListName, clusterVector);
        if (!clusterVector.empty())
        {
            HitType view{LArClusterHelper::GetClusterHitType(clusterVector.front())};
            clusterMap[view] = clusterVector;
        }
    }

    const PfoList *pTemporaryList{nullptr};
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTemporaryList, temporaryListName));

    bool clustersMatched{false};
    do
    {
        clustersMatched = false;
        float bestFom{std::numeric_limits<float>::max()};
        const Cluster *pBestClusterU{nullptr}, *pBestClusterV{nullptr}, *pBestClusterW{nullptr};
        for (const Cluster *const pClusterU : clusterMap[TPC_VIEW_U])
        {
            const TwoDSlidingFitResult fitResultU(pClusterU, m_fitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U));
            // Complex trajectories can have multiple fit segments - we're going to ignore those here
            if (fitResultU.GetFitSegmentList().size() != 1)
                continue;
            for (const Cluster *const pClusterV : clusterMap[TPC_VIEW_V])
            {
                const TwoDSlidingFitResult fitResultV(pClusterV, m_fitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V));
                if (fitResultV.GetFitSegmentList().size() != 1)
                    continue;

                for (const Cluster *const pClusterW : clusterMap[TPC_VIEW_W])
                {
                    const TwoDSlidingFitResult fitResultW(pClusterW, m_fitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W));
                    if (fitResultW.GetFitSegmentList().size() != 1)
                        continue;

                    if (!(PandoraContentApi::IsAvailable(*this, pClusterU) && PandoraContentApi::IsAvailable(*this, pClusterV) &&
                        PandoraContentApi::IsAvailable(*this, pClusterW)))
                        continue;
                    const float fom{this->GetFigureOfMerit(fitResultU, fitResultV, fitResultW)};
                    if (fom < bestFom)
                    {
                        bestFom = fom;
                        pBestClusterU = pClusterU;
                        pBestClusterV = pClusterV;
                        pBestClusterW = pClusterW;
                        clustersMatched = true;
                    }
                }
            }
        }
        
        if (clustersMatched)
        {
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            ClusterList clusterListU{{pBestClusterU}};
            ClusterList clusterListV{{pBestClusterV}};
            ClusterList clusterListW{{pBestClusterW}};
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListU, "U", RED));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListV, "V", GREEN));
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListW, "W", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

            PandoraContentApi::ParticleFlowObject::Parameters parameters;

            parameters.m_particleId = MU_MINUS;
            parameters.m_charge = PdgTable::GetParticleCharge(parameters.m_particleId.Get());
            parameters.m_mass = PdgTable::GetParticleMass(parameters.m_particleId.Get());
            parameters.m_energy = 0.f;
            parameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
            for (const Cluster *const pCluster : {pBestClusterU, pBestClusterV, pBestClusterW})
                parameters.m_clusterList.emplace_back(pCluster);
            const ParticleFlowObject *pPfo{nullptr};
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pPfo));
        }
    }
    while (clustersMatched);

    if (!pTemporaryList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*this, m_pfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, m_pfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float MyParticleCreationAlgorithm::GetFigureOfMerit(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV,
    const TwoDSlidingFitResult &fitResultW) const
{
    const FitSegment &fitSegmentU{fitResultU.GetFitSegmentList().front()};
    const FitSegment &fitSegmentV{fitResultV.GetFitSegmentList().front()};
    const FitSegment &fitSegmentW{fitResultW.GetFitSegmentList().front()};

    float chi2{0.f};
    const double xMin{std::max({fitSegmentU.GetMinX(), fitSegmentV.GetMinX(), fitSegmentW.GetMinX()})};
    const double xMax{std::min({fitSegmentU.GetMaxX(), fitSegmentV.GetMaxX(), fitSegmentW.GetMaxX()})};

    int consideredPoints{0};
    for (int n = 0; n < m_nSamplingPoints; ++n)
    {
        const double x{xMin + (xMax - xMin) * n / m_nSamplingPoints};
        CartesianVector fitVectorU(0, 0, 0), fitVectorV(0, 0, 0), fitVectorW(0, 0, 0);
        CartesianVector fitDirectionU(0, 0, 0), fitDirectionV(0, 0, 0), fitDirectionW(0, 0, 0);

        if (fitResultU.GetTransverseProjection(x, fitSegmentU, fitVectorU, fitDirectionU) != STATUS_CODE_SUCCESS ||
            fitResultV.GetTransverseProjection(x, fitSegmentV, fitVectorV, fitDirectionV) != STATUS_CODE_SUCCESS ||
            fitResultW.GetTransverseProjection(x, fitSegmentW, fitVectorW, fitDirectionW) != STATUS_CODE_SUCCESS)
            continue;

        const float u{fitVectorU.GetZ()}, v{fitVectorV.GetZ()}, w{fitVectorW.GetZ()};

        const float uv2w{LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, u, v)};
        const float uw2v{LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, u, w)};
        const float vw2u{LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, v, w)};

        const float deltaU{(vw2u - u) * fitDirectionU.GetX()};
        const float deltaV{(uw2v - v) * fitDirectionV.GetX()};
        const float deltaW{(uv2w - w) * fitDirectionW.GetX()};

        chi2 += deltaU * deltaU + deltaV * deltaV + deltaW * deltaW;
        ++consideredPoints;
    }

    return consideredPoints > 0 ? chi2 : std::numeric_limits<float>::max();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyParticleCreationAlgorithm::GetLongClusters(const std::string &clusterListName, ClusterVector &sortedClusters) const
{
    const ClusterList *pClusterList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() >= 3)
            sortedClusters.emplace_back(pCluster);
    }

    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FitWindow", m_fitWindow));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

