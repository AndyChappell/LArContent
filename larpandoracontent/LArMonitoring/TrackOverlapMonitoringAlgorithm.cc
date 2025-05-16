/**
 *  @file   larpandoracontent/LArMonitoring/TrackOverlapMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TrackOverlapMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

TrackOverlapMonitoringAlgorithm::TrackOverlapMonitoringAlgorithm() :
    m_visualise{true},
    m_writeFile{false},
    m_vertexRadius{10.f},
    m_distance{2.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapMonitoringAlgorithm::~TrackOverlapMonitoringAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootFileName, m_rootTreeName, "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::Run()
{
    if (m_visualise)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1, 1, 1));
    }

    this->CreateMaps();
    MCToMCMap overlapCandidates;
    this->FindTrueOverlapCandidates(overlapCandidates);


    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::CreateMaps()
{
    const PfoList *pfoList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pfoList));
    if (!pfoList || pfoList->empty())
        return STATUS_CODE_SUCCESS;

    for (const Pfo *const pPfo : *pfoList)
    {
        CaloHitList allCaloHits;
        LArPfoHelper::GetAllCaloHits(pPfo, allCaloHits);
        for (const CaloHit *const pCaloHit : allCaloHits)
        {
            const MCParticleWeightMap &contributionMap{pCaloHit->GetMCParticleWeightMap()};
            for (const auto &[pMC, weight] : contributionMap)
            {
                if (weight <= 0.f)
                    continue;

                m_pfoToMCMap[pPfo].insert(pMC);
                m_mcToPfoMap[pMC].insert(pPfo);
                m_mcToHitsMap[pMC].insert(pCaloHit);
            }
        }
    }

    for (const auto &[pMC, caloHits] : m_mcToHitsMap)
    {
        const int pdg{std::abs(pMC->GetParticleId())};
        if (caloHits.empty() || (pdg == 11|| pdg == 22) || pMC->GetMomentum().GetMagnitude() < 1e-3)
            continue;
        if (m_vertexToMCMap.find(pMC->GetVertex()) == m_vertexToMCMap.end())
            m_vertexToMCMap[pMC->GetVertex()] = MCParticleSet();
        m_vertexToMCMap[pMC->GetVertex()].insert(pMC);
        // Might want to think about endpoints too, but that makes direction comparison harder
        //m_vertexToMCMap[pMC->GetEndPoint()].insert(pMC);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::FindTrueOverlapCandidates(MCToMCMap &overlapCandidates) const
{
    const LArTransformationPlugin *pTransform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    for (const auto &[vertex, mcParticles] : m_vertexToMCMap)
    {
        if (mcParticles.size() < 2)
            continue;

        // Get all of the hits in each MC particle and find the closest approach for each hit in the other MC particles.
        // Retain hits within a couple of cm of a hit in each of the other MC particles.
        // If a sufficiently large fraction of hits are close within a 10cm radius of the vertex, we have an overlap candidate.
        for (auto iter1 = mcParticles.begin(); iter1 != mcParticles.end(); ++iter1)
        {
            const MCParticle *const pMC1{*iter1};
            CaloHitList uHits1, vHits1, wHits1;
            this->CollectHitsByView(pMC1, uHits1, vHits1, wHits1);
            float x{0}, u{0}, v{0}, w{0};
            (void)v;
            (void)w;
            LArVertexHelper::GetTrueVertexPosition(pMC1->GetVertex(), pTransform, x, u, v, w);
            Eigen::RowVector2f uVertex(x, u);
            Eigen::MatrixXf uHitMatrix1; // Size will be set by FilterHits
            this->VectorizeAndFilterHits(uHits1, uVertex, m_vertexRadius, uHitMatrix1);
            // Need to do the same filtering for the second batch of hits

            for (auto iter2 = std::next(iter1); iter2 != mcParticles.end(); ++iter2)
            {
                const MCParticle *const pMC2{*iter2};
                if (pMC1->GetMomentum().GetMagnitude() < 1e-3 || pMC2->GetMomentum().GetMagnitude() < 1e-3)
                    continue;

                // Veto particles with too large an opening angle - suggest being conservative here, mainly want to avoid back-to-back cases
                const CartesianVector dir1{pMC1->GetMomentum().GetUnitVector()};
                const CartesianVector dir2{pMC2->GetMomentum().GetUnitVector()};
                const float costheta{dir1.GetDotProduct(dir2)};
                if (costheta < 0.866f)
                    continue;

                CaloHitList uHits2, vHits2, wHits2;
                this->CollectHitsByView(pMC2, uHits2, vHits2, wHits2);
                Eigen::MatrixXf uHitMatrix2; // Size will be set by FilterHits
                this->VectorizeAndFilterHits(uHits2, uVertex, m_vertexRadius, uHitMatrix2);
                Eigen::MatrixXf uHitMatrix1Filtered, uHitMatrix2Filtered;
                this->GetDifferenceAndFilterHits(uHitMatrix1, uHitMatrix2, m_distance, uHitMatrix1Filtered, uHitMatrix2Filtered);

                // Find shared hits by hashing one lot of filtered hits and storing in a set, then look for hashed hits from the other lot of filtered
                // hits in the set
                std::unordered_set<MatrixHit, MatrixHitHash> uHitSet;
                for (int i = 0; i < uHitMatrix1Filtered.rows(); ++i)
                {
                    uHitSet.insert({uHitMatrix1Filtered(i, 0), uHitMatrix1Filtered(i, 1)});
                }
                size_t sharedHits{0};
                for (int i = 0; i < uHitMatrix2Filtered.rows(); ++i)
                {
                    if (uHitSet.find({uHitMatrix2Filtered(i, 0), uHitMatrix2Filtered(i, 1)}) != uHitSet.end())
                        ++sharedHits;
                }
                const size_t uniqueHits1{uHitMatrix1Filtered.rows() - sharedHits};
                const size_t uniqueHits2{uHitMatrix2Filtered.rows() - sharedHits};
                const bool isCandidate{sharedHits > 2 || (uniqueHits1 >= 2 && uniqueHits2 >= 2)};
                if (isCandidate)
                {
                    overlapCandidates[pMC1].insert(pMC2);
                    overlapCandidates[pMC2].insert(pMC1);
                }
                if (m_visualise && isCandidate && uHitMatrix1Filtered.rows() > 0 && uHitMatrix2Filtered.rows() > 0)
                {
                    const CartesianVector pos(x, 0, u); (void)pos;
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits1, "1 (" + std::to_string(uHitMatrix1.rows()) + ")", BLACK));
                    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits2, "2 (" + std::to_string(uHitMatrix2.rows()) + ")", RED));
                    for (int i = 0; i < uHitMatrix1Filtered.rows(); ++i)
                    {
                        const CartesianVector pos1{uHitMatrix1Filtered(i, 0), 0, uHitMatrix1Filtered(i, 1)}; (void)pos1;
                        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos1, "1", BLACK, 1));
                    }
                    for (int i = 0; i < uHitMatrix2Filtered.rows(); ++i)
                    {
                        const CartesianVector pos2{uHitMatrix2Filtered(i, 0), 0, uHitMatrix2Filtered(i, 1)}; (void)pos2;
                        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos2, "2", RED, 1));
                    }

                    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos, "vtx", BLUE, 1));
                    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                }
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::AssessPfos(const MCToMCMap &overlapCandidates, const PfoList &pfoList) const
{
    (void)overlapCandidates;
    (void)pfoList;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackOverlapMonitoringAlgorithm::CollectHitsByView(const MCParticle *const pMC, CaloHitList &uHits, CaloHitList &vHits, CaloHitList &wHits) const
{
    const CaloHitSet &caloHits{m_mcToHitsMap.at(pMC)};
    for (const CaloHit *const pCaloHit : caloHits)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                uHits.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_V:
                vHits.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_W:
                wHits.emplace_back(pCaloHit);
                break;
            default:
                break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackOverlapMonitoringAlgorithm::VectorizeAndFilterHits(const pandora::CaloHitList &hits, const Eigen::RowVector2f &vertex, const float distance, Eigen::MatrixXf &filteredHits) const
{
    const float threshold{distance * distance};
    Eigen::MatrixXf allHitsMatrix(hits.size(), 2);
    LArEigenHelper::Vectorize(hits, allHitsMatrix);

    Eigen::MatrixXf diffsSq{(allHitsMatrix.rowwise() - vertex).rowwise().squaredNorm()};
    std::vector<int> keepIndices;
    for (int i = 0; i < diffsSq.rows(); ++i)
    {
        if (diffsSq(i, 0) < threshold)
            keepIndices.emplace_back(i);
    }
    filteredHits.resize(keepIndices.size(), allHitsMatrix.cols());
    for (int i = 0; i < static_cast<int>(keepIndices.size()); ++i)
        filteredHits.row(i) = allHitsMatrix.row(keepIndices[i]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackOverlapMonitoringAlgorithm::GetDifferenceAndFilterHits(const Eigen::MatrixXf &hits1, const Eigen::MatrixXf &hits2, const float distance, Eigen::MatrixXf &filteredHits1, Eigen::MatrixXf &filteredHits2) const
{
    const float threshold{distance * distance};
    Eigen::VectorXf sqNorm1(hits1.rowwise().squaredNorm());
    Eigen::VectorXf sqNorm2(hits2.rowwise().squaredNorm());
    Eigen::MatrixXf distances = sqNorm1.replicate(1, hits2.rows()) + sqNorm2.transpose().replicate(hits1.rows(), 1) - 2.f * (hits1 * hits2.transpose());

    std::unordered_set<int> keepIndices1, keepIndices2;
    for (int i = 0; i < distances.rows(); ++i)
    {
        for (int j = 0; j < distances.cols(); ++j)
        {
            if (distances(i, j) < threshold)
            {
                keepIndices1.insert(i);
                keepIndices2.insert(j);
            }
        }
    }

    filteredHits1.resize(keepIndices1.size(), hits1.cols());
    auto iter1 = keepIndices1.begin();
    for (int i = 0; i < static_cast<int>(keepIndices1.size()); ++i, ++iter1)
        filteredHits1.row(i) = hits1.row(*iter1);
    filteredHits2.resize(keepIndices2.size(), hits2.cols());
    auto iter2 = keepIndices2.begin();
    for (int i = 0; i < static_cast<int>(keepIndices2.size()); ++i, ++iter2)
        filteredHits2.row(i) = hits2.row(*iter2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackOverlapMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));

    if (m_writeFile)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRadius", m_vertexRadius));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OverlapDistance", m_distance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
