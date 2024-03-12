/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlClusterAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning clustering algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArTwoDReco/DlClusterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArDelaunayTriangulationHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include <random>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

StatusCode DlClusterAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (pCaloHitList->empty())
            continue;
        std::cout << "Num hits: " << pCaloHitList->size() << std::endl;
/*        CaloHitList subsetList;
        auto iter{pCaloHitList->begin()};
        subsetList.emplace_back(*iter);
        std::advance(iter, 2);
        subsetList.emplace_back(*iter);
        std::advance(iter, 10);
        subsetList.emplace_back(*iter);
        std::advance(iter, 20);
        subsetList.emplace_back(*iter);*/

        LArDelaunayTriangulationHelper::VertexVector vertices;
        LArDelaunayTriangulationHelper::TriangleVector triangles;
        LArDelaunayTriangulationHelper::Triangulate(*pCaloHitList, vertices, triangles);
        std::map<const LArDelaunayTriangulationHelper::Vertex*, LArDelaunayTriangulationHelper::EdgeVector> vertexEdgeMap;

        LArDelaunayTriangulationHelper::EdgeVector edges;
        std::random_device device;
        std::mt19937 generator(device());
        std::bernoulli_distribution coinFlip(1);
        for (const LArDelaunayTriangulationHelper::Triangle *pTriangle : triangles)
        {
            if (coinFlip(generator))
            {
                LArDelaunayTriangulationHelper::EdgeVector triEdges;
                pTriangle->GetEdges(triEdges);
                const float len0Square{triEdges.at(0)->m_v0->DistanceSquared(*triEdges.at(0)->m_v1)};
                const float len1Square{triEdges.at(1)->m_v0->DistanceSquared(*triEdges.at(1)->m_v1)};
                const float len2Square{triEdges.at(2)->m_v0->DistanceSquared(*triEdges.at(2)->m_v1)};
                const float minLen{std::min({len0Square, len1Square, len2Square})};
                const float maxLen{std::max({len0Square, len1Square, len2Square})};
                if ((m_maxEdgeRatioSquared * minLen) < maxLen)
                {
                    // Only base of triangle is a reasonable link - delete the other edges
                    if (len2Square > minLen || len2Square > m_maxEdgeLengthSquared)
                        triEdges.erase(triEdges.begin() + 2);
                    if (len1Square > minLen || len1Square > m_maxEdgeLengthSquared)
                        triEdges.erase(triEdges.begin() + 1);
                    if (len0Square > minLen || len0Square > m_maxEdgeLengthSquared)
                        triEdges.erase(triEdges.begin() + 0);
                }
                else
                {
                    // Remove excessively long edges
                    if (len2Square > m_maxEdgeLengthSquared)
                        triEdges.erase(triEdges.begin() + 2);
                    if (len1Square > m_maxEdgeLengthSquared)
                        triEdges.erase(triEdges.begin() + 1);
                    if (len0Square > m_maxEdgeLengthSquared)
                        triEdges.erase(triEdges.begin() + 0);
                }
                for (const LArDelaunayTriangulationHelper::Edge *pEdge : triEdges)
                {
                    if (std::find_if(edges.begin(), edges.end(),
                        [pEdge](const LArDelaunayTriangulationHelper::Edge *const &edge) { return *edge == *pEdge; }) == edges.end())
                    {
                        edges.emplace_back(pEdge);
                    }
                }
            }
        }
        if (m_prune)
        {
            LArDelaunayTriangulationHelper::EdgeVector prunedEdges;
            for (const LArDelaunayTriangulationHelper::Edge *pEdge : edges)
            {
                vertexEdgeMap[pEdge->m_v0].emplace_back(pEdge);
                vertexEdgeMap[pEdge->m_v1].emplace_back(pEdge);
            }
            for (const LArDelaunayTriangulationHelper::Edge *pEdge1 : edges)
            {
                const float thisLength{pEdge1->LengthSquared()};
                bool shortest{true};
                for (const LArDelaunayTriangulationHelper::Edge *pEdge2 : vertexEdgeMap[pEdge1->m_v0])
                {
                    if (pEdge1 == pEdge2)
                        continue;
                    const float compareLength{pEdge2->LengthSquared()};
                    if (compareLength < thisLength)
                    {
                        shortest = false;
                        break;
                    }
                }
                if (!shortest)
                {
                    shortest = true;
                    for (const LArDelaunayTriangulationHelper::Edge *pEdge2 : vertexEdgeMap[pEdge1->m_v1])
                    {
                        if (pEdge1 == pEdge2)
                            continue;
                        const float compareLength{pEdge2->LengthSquared()};
                        if (compareLength < thisLength)
                        {
                            shortest = false;
                            break;
                        }
                    }
                }
                if (shortest)
                    prunedEdges.emplace_back(pEdge1);
            }
            edges = prunedEdges;
        }
        for (const LArDelaunayTriangulationHelper::Edge *pEdge : edges)
        {
            CartesianVector start(pEdge->m_v0->m_x, 0, pEdge->m_v0->m_z);
            CartesianVector end(pEdge->m_v1->m_x, 0, pEdge->m_v1->m_z);
            if (pEdge->m_v0->m_pCaloHit && pEdge->m_v1->m_pCaloHit)
            {
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &end, "e", BLUE, 1, 1));
            }
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        // Visualise the triangulation
        // COnvert to graph
        // Visualise the graph

        for (const LArDelaunayTriangulationHelper::Vertex *pVertex : vertices)
            delete pVertex;
        for (const LArDelaunayTriangulationHelper::Triangle *pTriangle : triangles)
            delete pTriangle;
    }

    return STATUS_CODE_SUCCESS;
    //return this->PrepareTrainingSample();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterAlgorithm::PrepareTrainingSample()
{
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (pCaloHitList->empty())
            continue;
        CaloHitVector caloHitVector(pCaloHitList->begin(), pCaloHitList->end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        std::map<const MCParticle *, CaloHitList> mcToHitsMap;
        for (const CaloHit *pCaloHit : caloHitVector)
        {
            try
            {
                const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                mcToHitsMap[pMC].emplace_back(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
            }
        }

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(mcToHitsMap.size()));
        for (const auto & [pMC, caloHitList] : mcToHitsMap)
        {
            (void)pMC;
            featureVector.emplace_back(static_cast<double>(caloHitList.size()));
            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                const CartesianVector &position{pCaloHit->GetPositionVector()};
                const float x{position.GetX()}, z{position.GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};
                featureVector.emplace_back(x);
                featureVector.emplace_back(z);
                featureVector.emplace_back(adc);
            }
        }
        const HitType view{pCaloHitList->front()->GetHitType()};
        std::string suffix{view == TPC_VIEW_U ? "_U.csv" : view == TPC_VIEW_V ? "_V.csv" : "_W.csv"};
        LArMvaHelper::ProduceTrainingExample(m_outputFilePrefix + suffix, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFilePrefix", m_outputFilePrefix));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MaxEdgeLength", m_maxEdgeLengthSquared));
    m_maxEdgeLengthSquared *= m_maxEdgeLengthSquared;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MaxEdgeRatio", m_maxEdgeRatioSquared));
    m_maxEdgeRatioSquared *= m_maxEdgeRatioSquared;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Prune", m_prune));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
