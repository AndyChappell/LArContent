/**
 *  @file   larpandoracontent/LArHelpers/LArTopologyHelper.cc
 *
 *  @brief  Implementation of the principal curve analysis helper class.
 *
 *  $Log: $
 */

#include "Pandora/PandoraEnumeratedTypes.h"
#include "Pandora/PdgTable.h"

#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArTopologyHelper.h"

#include <Eigen/Dense>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

LArTopologyHelper::Filter::Filter(const size_t minHits, const bool neutrons, const bool deltas, const bool compton, const bool photoElectric) :
    m_minHits(minHits),
    m_neutrons(neutrons),
    m_deltas(deltas),
    m_compton(compton),
    m_photoElectric(photoElectric)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::GetTopologicalVertices(const CaloHitList &caloHitList, const Filter &filter,
    const LArTransformationPlugin *const pTransform, CartesianPointVector &vertices)
{
    if (caloHitList.empty())
        return;
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    LArMCParticleHelper::GetMCToHitsMap(caloHitList, mcToHitsMap, true);
    LArTopologyHelper::FilterVertices(mcToHitsMap, filter);
    MCVertexMap mcToVertexMap;
    for (const auto &[pMC, mcHits] : mcToHitsMap)
    {
        const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (!pLArMC)
            continue;
        std::cout << "MC: " << pMC->GetParticleId() << " (" << mcHits.size() << " hits)" << " Process: " << pLArMC->GetProcess() << std::endl;
        LArTopologyHelper::GetProjectedTrueVertices(pMC, caloHitList.front()->GetHitType(), pTransform, mcToVertexMap);
    }

    if (!mcToVertexMap.empty())
    {
        // Find the closest hit to each vertex, along with the distance between the hit and the vertex
        MCVertexMap mcToMatchedVertexMap;
        LArTopologyHelper::MatchHitToVertex(mcToHitsMap, mcToVertexMap, mcToMatchedVertexMap);

        // Now consolidate any duplicate vertices and associate the relevant hits to the consolidated vertex
        VertexHitsMap vertexHitsMap;
        LArTopologyHelper::ConsolidateVertices(mcToMatchedVertexMap, mcToHitsMap, vertexHitsMap);

        for (const auto &[vertex, _] : vertexHitsMap)
            vertices.emplace_back(vertex);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::GetProjectedTrueVertices(const MCParticle *const pMC, const HitType view, const LArTransformationPlugin *const pTransform,
    MCVertexMap &mcVertexMap)
{
    CartesianPointVector mcVertices;
    if (pMC->GetParticleId() == PHOTON)
    {
        // Photon hit deposition begins at the endpoint (i.e. after conversion)
        mcVertices.emplace_back(pMC->GetEndpoint());
    }
    else
    {
        // For track-like particles use both the start and end points
        mcVertices.emplace_back(pMC->GetVertex());
        if (std::abs(pMC->GetParticleId()) != E_MINUS)
        {
            mcVertices.emplace_back(pMC->GetEndpoint());
        }
        else
        {
            const LArMCParticle *const pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
            if (pLArMC && (pLArMC->GetProcess() == MC_PROC_DECAY))
                mcVertices.emplace_back(pMC->GetEndpoint());
        }
    }
    for (const CartesianVector &mcVertex : mcVertices)
    {
        CartesianVector viewVertex(0, 0, 0);
        switch (view)
        {
            case TPC_VIEW_U:
                viewVertex.SetValues(mcVertex.GetX(), 0, pTransform->YZtoU(mcVertex.GetY(), mcVertex.GetZ()));
                break;
            case TPC_VIEW_V:
                viewVertex.SetValues(mcVertex.GetX(), 0, pTransform->YZtoV(mcVertex.GetY(), mcVertex.GetZ()));
                break;
            case TPC_VIEW_W:
                viewVertex.SetValues(mcVertex.GetX(), 0, pTransform->YZtoW(mcVertex.GetY(), mcVertex.GetZ()));
                break;
            default:
                break;
        }
        mcVertexMap[pMC].emplace_back(viewVertex);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::MatchHitToVertex(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, const MCVertexMap &mcVertexMap,
    MCVertexMap &mcToMatchedVertexMap)
{
    for (const auto &[pMC, caloHitList] : mcToHitsMap)
    {
        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        for (const CartesianVector &vertex : mcVertexMap.at(pMC))
        {
            const CartesianPointVector vertices(1, vertex);
            Eigen::MatrixXf hitMat(caloHitList.size(), 2), vertMat(1, 2);
            LArEigenHelper::Vectorize(caloHitList, hitMat);
            LArEigenHelper::Vectorize({vertices}, vertMat);

            if (hitMat.cols() != vertMat.cols())
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            Eigen::VectorXf hitSq{hitMat.rowwise().squaredNorm()};
            Eigen::VectorXf vertSq{vertMat.rowwise().squaredNorm()};

            // Compute all pairwise squared distances:
            // D2(i,j) = |h_i|^2 + |v_j|^2 - 2 * h_i.dot(v_j)
            Eigen::MatrixXf dot{hitMat * vertMat.transpose()};
            Eigen::MatrixXf d2{(-2.0f * dot).colwise() + hitSq}; // add hit norms per column
            d2 = d2.rowwise() + vertSq.transpose();              // add vertex norms per row

            Eigen::VectorXf col{d2.col(0)};
            Eigen::Index idx{0};
            col.minCoeff(&idx);

            const int index{static_cast<int>(idx)};
            mcToMatchedVertexMap[pMC].emplace_back(caloHitVector[index]->GetPositionVector());
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::ConsolidateVertices(const MCVertexMap &mcToMatchedVertexMap, const LArMCParticleHelper::MCContributionMap &mcToHitsMap, VertexHitsMap &vertexHitsMap)
{
    // We only want the leading vertex in an EM cascade
    MCVertexMap filteredMCToMatchedVertexMap;
    for (const auto &[pMC, vertices] : mcToMatchedVertexMap)
    {
        const int pdg{std::abs(pMC->GetParticleId())};
        if (pdg == E_MINUS)
        {
            const MCParticle *pParent{pMC};
            bool foundElectronParent{false};
            while (!pParent->GetParentList().empty())
            {
                pParent = pParent->GetParentList().front();
                if (std::abs(pParent->GetParticleId()) == E_MINUS)
                {
                    foundElectronParent = true;
                    break;
                }
            }
            if (!foundElectronParent)
                filteredMCToMatchedVertexMap[pMC] = vertices;
        }
        else
        {
            filteredMCToMatchedVertexMap[pMC] = vertices;
        }
    }

    for (const auto &[pMC, vertices] : filteredMCToMatchedVertexMap)
    {
        bool skipEndpoint{false};
        const MCParticleList &children{pMC->GetDaughterList()};
        // Look for children that have already been matched to a vertex - if so, skip this endpoint
        for (const MCParticle *const pChild : children)
        {
            // We still want the endpoint if the child is a photon, and therefore offset from the track end
            if (std::abs(pChild->GetParticleId()) == PHOTON)
                continue;
            if (mcToMatchedVertexMap.find(pChild) != mcToMatchedVertexMap.end())
            {
                skipEndpoint = true;
                break;
            }
        }
        for (const CartesianVector &vertex : vertices)
        {
            if (vertexHitsMap.find(vertex) == vertexHitsMap.end())
                vertexHitsMap[vertex] = CaloHitList();
            const CaloHitList &caloHitList{mcToHitsMap.at(pMC)};
            vertexHitsMap[vertex].insert(vertexHitsMap[vertex].end(), caloHitList.begin(), caloHitList.end());
            if (skipEndpoint)
                break;
        }
    }

    // Finally, if multiple vetices are found within a small distance, pick the vertex from the highest energy particle
    const float minDistanceSquared{2.f * 2.f}; // 2 cm
    CartesianPointVector verticesToRemove;
    for (auto it1 = vertexHitsMap.begin(); it1 != vertexHitsMap.end(); ++it1)
    {
        const CartesianVector &vertex1{it1->first};
        for (auto it2 = std::next(it1); it2 != vertexHitsMap.end(); ++it2)
        {
            const CartesianVector &vertex2{it2->first};
            const float distanceSquared{(vertex1 - vertex2).GetMagnitudeSquared()};
            if (distanceSquared < minDistanceSquared)
            {
                // Find the total energy for each vertex
                float energy1{0.f}, energy2{0.f};
                for (const CaloHit *const pCaloHit : it1->second)
                    energy1 += pCaloHit->GetInputEnergy();
                for (const CaloHit *const pCaloHit : it2->second)
                    energy2 += pCaloHit->GetInputEnergy();
                // Erase the lower energy vertex - if they're the same, likely vertex and endpoint from one particle, keep both
                if (energy1 > energy2)
                    verticesToRemove.emplace_back(vertex2);
                else if (energy2 > energy1)
                    verticesToRemove.emplace_back(vertex1);
            }
        }
    }
    for (const CartesianVector &vertex : verticesToRemove)
        vertexHitsMap.erase(vertex);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::FilterVertices(LArMCParticleHelper::MCContributionMap &mcToHitsMap, const Filter &filter)
{
    if (filter.m_minHits > 0)
        LArTopologyHelper::FilterByMinHits(mcToHitsMap, filter.m_minHits);
    if (filter.m_deltas)
        LArTopologyHelper::FilterDeltaRays(mcToHitsMap);
    if (filter.m_neutrons)
        LArTopologyHelper::FilterNeutronInducedParticles(mcToHitsMap);
    if (filter.m_compton)
        LArTopologyHelper::FilterComptonScattering(mcToHitsMap);
    if (filter.m_photoElectric)
        LArTopologyHelper::FilterPhotoElectrons(mcToHitsMap);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::FilterDeltaRays(LArMCParticleHelper::MCContributionMap &mcToHitsMap)
{
    for (auto it = mcToHitsMap.begin(); it != mcToHitsMap.end();)
    {
        const MCParticle *const pMC{it->first};
        if (std::abs(pMC->GetParticleId()) != E_MINUS)
        {
            ++it;
            continue;
        }
        // Found electron, now check if it is a descendant of an ionisation process
        const MCParticle *pCurrentMC{pMC};
        bool foundDelta{false};
        while (!pCurrentMC->GetParentList().empty())
        {
            if (LArMCParticleHelper::IsIonisation(pCurrentMC))
            {
                it = mcToHitsMap.erase(it);
                foundDelta = true;
                break;
            }
            pCurrentMC = pCurrentMC->GetParentList().front();
        }
        if (!foundDelta)
            ++it;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::FilterNeutronInducedParticles(LArMCParticleHelper::MCContributionMap &mcToHitsMap)
{
    for (auto it = mcToHitsMap.begin(); it != mcToHitsMap.end();)
    {
        const MCParticle *const pMC{it->first};
        bool foundNeutronParent{false};
        if (std::abs(pMC->GetParticleId()) == PHOTON)
        {
            const MCParticle *pCurrentMC{pMC};
            while (!pCurrentMC->GetParentList().empty())
            {
                pCurrentMC = pCurrentMC->GetParentList().front();
                if (std::abs(pCurrentMC->GetParticleId()) == NEUTRON)
                {
                    it = mcToHitsMap.erase(it);
                    foundNeutronParent = true;
                    break;
                }
            }
        }
        if (!foundNeutronParent)
            ++it;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::FilterComptonScattering(LArMCParticleHelper::MCContributionMap &mcToHitsMap)
{
    for (auto it = mcToHitsMap.begin(); it != mcToHitsMap.end();)
    {
        const MCParticle *const pMC{it->first};
        if (std::abs(pMC->GetParticleId()) != E_MINUS)
        {
            ++it;
            continue;
        }
        // Found electron, now check for Compton scatter
        const LArMCParticle *pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC && pLArMC->GetProcess() == MC_PROC_COMPT)
            it = mcToHitsMap.erase(it);
        else
            ++it;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::FilterPhotoElectrons(LArMCParticleHelper::MCContributionMap &mcToHitsMap)
{
    for (auto it = mcToHitsMap.begin(); it != mcToHitsMap.end();)
    {
        const MCParticle *const pMC{it->first};
        if (std::abs(pMC->GetParticleId()) != E_MINUS)
        {
            ++it;
            continue;
        }
        // Found electron, now check for photoelectric effect
        const LArMCParticle *pLArMC{dynamic_cast<const LArMCParticle *>(pMC)};
        if (pLArMC && pLArMC->GetProcess() == MC_PROC_PHOT)
            it = mcToHitsMap.erase(it);
        else
            ++it;
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArTopologyHelper::FilterByMinHits(LArMCParticleHelper::MCContributionMap &mcToHitsMap, const size_t minHits)
{
    for (auto it = mcToHitsMap.begin(); it != mcToHitsMap.end();)
    {
        const MCParticle *const pMC{it->first};
        const int pdg{std::abs(pMC->GetParticleId())};
        // Don't filter photons or electrons - this is handled by other, more appropriate mechanisms
        if (pdg == PHOTON || pdg == E_MINUS)
        {
            ++it;
            continue;
        }
        const CaloHitList &caloHitList{it->second};
        if (caloHitList.size() < minHits)
            it = mcToHitsMap.erase(it);
        else
            ++it;
    }
}

} // namespace lar_content
