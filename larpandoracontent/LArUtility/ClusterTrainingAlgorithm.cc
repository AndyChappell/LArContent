/**
 *  @file   larpandoracontent/LArUtility/ClusterTrainingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster training algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArTpcGeometryHelper.h"
#include "larpandoracontent/LArUtility/ClusterTrainingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterTrainingAlgorithm::ClusterTrainingAlgorithm() : m_initialisation{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterTrainingAlgorithm::Run()
{
    const LArTransformationPlugin *const pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    LArTpcGeometryHelper helper{LArTpcGeometryHelper::GetInstance(this, pTransform)};
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));
    for (unsigned int t = 0; t < 2; ++t)
    {
        for (unsigned int c = t % 2; c < 24; c += 2)
        {
            LArTpcGeometryHelper::VolumeId id(0, t, c);
            TpcHitVolume &volume{helper.GetTpcHitVolume(id)};
            if (m_initialisation)
            {
                volume.Add(*pCaloHitList);
                CaloHitList volumeHitsU, volumeHitsV, volumeHitsW;
                volume.GetHitList(TPC_VIEW_U, volumeHitsU);
                volume.GetHitList(TPC_VIEW_V, volumeHitsV);
                volume.GetHitList(TPC_VIEW_W, volumeHitsW);
                this->CreateHitList(volumeHitsU, "CaloHitListU_0_" + std::to_string(t) + "_" + std::to_string(c));
                this->CreateHitList(volumeHitsV, "CaloHitListV_0_" + std::to_string(t) + "_" + std::to_string(c));
                this->CreateHitList(volumeHitsW, "CaloHitListW_0_" + std::to_string(t) + "_" + std::to_string(c));
            }
            else
            {
                // Not sure why I need to run this line again, wouldn't have expected to, as singleton seems to persist as expected
                volume.Add(*pCaloHitList);
                TpcHitVolume::ClusterMap signalMap, backgroundMap;
                volume.GetCombinatorics(signalMap, backgroundMap);
                for (const auto &[ pCluster1, clusterList ] : signalMap)
                {
                    HitType view1{LArClusterHelper::GetClusterHitType(pCluster1)};
                    float xMin1{0.f}, xMax1{0.f};
                    pCluster1->GetClusterSpanX(xMin1, xMax1);
                    for (const Cluster *pCluster2 : clusterList)
                    {
                        float xMin2{0.f}, xMax2{0.f};
                        pCluster2->GetClusterSpanX(xMin2, xMax2);
                        const float xMin{std::max(xMin1, xMin2)}, xMax{std::min(xMax1, xMax2)};

                        if (xMin <= xMax)
                        {
                            HitType view2{LArClusterHelper::GetClusterHitType(pCluster2)};
                            CartesianPointVector hits1;
                            volume.GetLocalCoordinates(pCluster1, xMin, xMax, hits1);
                            CartesianPointVector hits2;
                            volume.GetLocalCoordinates(pCluster2, xMin, xMax, hits2);
                            if (hits1.size() >= 5 && hits2.size() >= 5)
                                this->WriteFeatureVector(view1, hits1, view2, hits2, true);
                        }
                    }
                }

                for (const auto &[ pCluster1, clusterList ] : backgroundMap)
                {
                    HitType view1{LArClusterHelper::GetClusterHitType(pCluster1)};
                    float xMin1{0.f}, xMax1{0.f};
                    pCluster1->GetClusterSpanX(xMin1, xMax1);
                    for (const Cluster *pCluster2 : clusterList)
                    {
                        float xMin2{0.f}, xMax2{0.f};
                        pCluster2->GetClusterSpanX(xMin2, xMax2);
                        const float xMin{std::max(xMin1, xMin2)}, xMax{std::min(xMax1, xMax2)};

                        if (xMin <= xMax)
                        {
                            HitType view2{LArClusterHelper::GetClusterHitType(pCluster2)};
                            CartesianPointVector hits1;
                            volume.GetLocalCoordinates(pCluster1, xMin, xMax, hits1);
                            CartesianPointVector hits2;
                            volume.GetLocalCoordinates(pCluster2, xMin, xMax, hits2);
                            if (hits1.size() >= 5 && hits2.size() >= 5)
                                this->WriteFeatureVector(view1, hits1, view2, hits2, false);
                        }
                    }
                }
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterTrainingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Initialise", m_initialisation));
    if (!m_initialisation)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingFilename", m_trainingFilename));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterTrainingAlgorithm::CreateHitList(const CaloHitList &caloHitList, std::string listName) const
{
    return PandoraContentApi::SaveList<CaloHitList>(*this, caloHitList, listName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterTrainingAlgorithm::WriteFeatureVector(const HitType view1, const CartesianPointVector &points1, const HitType view2,
    const CartesianPointVector &points2, const bool isSignal) const
{
    LArMvaHelper::MvaFeatureVector featureVector;

    featureVector.emplace_back(static_cast<double>(view1));
    featureVector.emplace_back(static_cast<double>(points1.size()));
    for (const CartesianVector &vec : points1)
    {
        featureVector.emplace_back(static_cast<double>(vec.GetX()));
        featureVector.emplace_back(static_cast<double>(vec.GetZ()));
    }
    featureVector.emplace_back(static_cast<double>(view2));
    featureVector.emplace_back(static_cast<double>(points2.size()));
    for (const CartesianVector &vec : points2)
    {
        featureVector.emplace_back(static_cast<double>(vec.GetX()));
        featureVector.emplace_back(static_cast<double>(vec.GetZ()));
    }

    return LArMvaHelper::ProduceTrainingExample(m_trainingFilename, isSignal, featureVector);
}

} // namespace lar_content

