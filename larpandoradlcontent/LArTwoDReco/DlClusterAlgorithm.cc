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

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

StatusCode DlClusterAlgorithm::Run()
{
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
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            vertices.emplace_back(new LArDelaunayTriangulationHelper::Vertex(pCaloHit));
        }

        LArDelaunayTriangulationHelper::TriangleVector triangles;
        const LArDelaunayTriangulationHelper::Triangle *bounds{LArDelaunayTriangulationHelper::MakeInitialBoundingTriangle(vertices)};
        triangles.emplace_back(bounds);
        for (const LArDelaunayTriangulationHelper::Vertex *const pVertex : vertices)
        {
            LArDelaunayTriangulationHelper::AddVertex(pVertex, triangles);
        }
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

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
