/**
 *  @file   larpandoracontent/LArHelpers/LArDelaunayTriangulationHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArDelaunayTriangulationHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

LArDelaunayTriangulationHelper::Triangle *LArDelaunayTriangulationHelper::MakeInitialBoundingTriangle(const VertexVector &vertices)
{
    float xMin{std::numeric_limits<float>::max()}, xMax{-std::numeric_limits<float>::max()}, zMin{std::numeric_limits<float>::max()},
        zMax{-std::numeric_limits<float>::max()};

    for (const Vertex *const pVertex : vertices)
    {
        xMin = std::min(xMin, pVertex->m_x);
        xMax = std::max(xMax, pVertex->m_x);
        zMin = std::min(zMin, pVertex->m_z);
        zMax = std::max(zMax, pVertex->m_z);
    }
    xMin -= 1;
    xMax += 1;
    zMin -= 1;
    zMax += 1;
    const float dx{xMax - xMin}, dz{zMax - zMin};
    const Vertex *const v0{new Vertex(0.5f * (xMin + xMax), zMin - dz)};
    const Vertex *const v1{new Vertex(v0->m_x - 1.5f * dx, v0->m_z + 3.f * dz)};
    const Vertex *const v2{new Vertex(v0->m_x + 1.5f * dx, v0->m_z + 3.f * dz)};

    return new Triangle(v0, v1, v2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDelaunayTriangulationHelper::AddVertex(const Vertex *const pVertex, TriangleVector &triangles)
{
    EdgeVector edges;
    // Remove triangles whose circumcircles contain the vertex and collect the unique edges
    for (auto iter = triangles.begin(); iter != triangles.end(); )
    {
        const Triangle *const triangle{*iter};
        if (triangle->IsInCircumCircle(pVertex))
        {
            triangle->GetUniqueEdges(edges);
            iter = triangles.erase(iter);
            std::cout << "Delete triangle: ";
            std::cout << "(" << triangle->m_v0->m_x << " " << triangle->m_v0->m_z << ") " <<
                "(" << triangle->m_v1->m_x << " " << triangle->m_v1->m_z << ") " <<
                "(" << triangle->m_v2->m_x << " " << triangle->m_v2->m_z << ") " << std::endl;
            delete triangle;
        }
        else
        {
            ++iter;
            std::cout << "Keep triangle: ";
            std::cout << "(" << triangle->m_v0->m_x << " " << triangle->m_v0->m_z << ") " <<
                "(" << triangle->m_v1->m_x << " " << triangle->m_v1->m_z << ") " <<
                "(" << triangle->m_v2->m_x << " " << triangle->m_v2->m_z << ") " << std::endl;
        }
    }

    for (const Edge *const pEdge : edges)
    {
        std::cout << "Edge: (" << pEdge->m_v0->m_x << " " << pEdge->m_v0->m_z << ") (" << pEdge->m_v1->m_x << " " << pEdge->m_v1->m_z << ")" << std::endl;
    }

    // Construct new triangles from the vertex and the unique edges
    for (const Edge *const pEdge : edges)
    {
        const Triangle *const triangle{new Triangle(pEdge->m_v0, pEdge->m_v1, pVertex)};
        std::cout << "New triangle: ";
        std::cout << "(" << triangle->m_v0->m_x << " " << triangle->m_v0->m_z << ") " <<
            "(" << triangle->m_v1->m_x << " " << triangle->m_v1->m_z << ") " <<
            "(" << triangle->m_v2->m_x << " " << triangle->m_v2->m_z << ") " << std::endl;

        triangles.emplace_back(triangle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDelaunayTriangulationHelper::ShrinkWrap(const Triangle &bounds, TriangleVector &triangles)
{
    EdgeVector boundingEdges;
    bounds.GetEdges(boundingEdges);
    for (auto iter = triangles.begin(); iter != triangles.end(); )
    {
        const Triangle *const triangle{*iter};
        if (triangle->SharesEdge(bounds))
        {
            iter = triangles.erase(iter);
            delete triangle;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Vertex::Vertex(const CaloHit *const pCaloHit) :
    m_x{pCaloHit->GetPositionVector().GetX()},
    m_z{pCaloHit->GetPositionVector().GetZ()},
    m_pCaloHit{pCaloHit}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Vertex::Vertex(const float x, const float z) :
    m_x{x},
    m_z{z},
    m_pCaloHit{nullptr}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Edge::Edge(const Vertex *const pVertex0, const Vertex *const pVertex1) :
    m_v0{pVertex0},
    m_v1{pVertex1}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Triangle::Triangle(const Vertex *const pVertex0, const Vertex *const pVertex1, const Vertex *const pVertex2) :
    m_v0{pVertex0},
    m_v1{pVertex1},
    m_v2{pVertex2},
    m_e01{new Edge(this->m_v0, this->m_v1)},
    m_e12{new Edge(this->m_v1, this->m_v2)},
    m_e20{new Edge(this->m_v2, this->m_v0)},
    m_ccx{0.f},
    m_ccz{0.f},
    m_ccr2{0.f}
{
    this->CalculateCircumCircle();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDelaunayTriangulationHelper::Triangle::IsInCircumCircle(const Vertex *const pVertex) const
{
    const float dx{pVertex->m_x - this->m_ccx};
    const float dz{pVertex->m_z - this->m_ccz};

    return (dx * dx + dz * dz) <= this->m_ccr2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDelaunayTriangulationHelper::Triangle::CalculateCircumCircle()
{
    const float minx{std::min({this->m_v0->m_x, this->m_v1->m_x, this->m_v2->m_x})};
    const float minz{std::min({this->m_v0->m_z, this->m_v1->m_z, this->m_v2->m_z})};
    const float maxx{std::max({this->m_v0->m_x, this->m_v1->m_x, this->m_v2->m_x})};
    const float maxz{std::max({this->m_v0->m_z, this->m_v1->m_z, this->m_v2->m_z})};
    const float ox{0.5f * (minx + maxx)};
    const float oz{0.5f * (minz + maxz)};

    const float dx0{this->m_v0->m_x - ox};
    const float dz0{this->m_v0->m_z - oz};
    const float dx1{this->m_v1->m_x - ox};
    const float dz1{this->m_v1->m_z - oz};
    const float dx2{this->m_v2->m_x - ox};
    const float dz2{this->m_v2->m_z - oz};
    const float scale{2 * (dx0 * (dz1 - dz2) + dx1 * (dz2 - dz0) + dx2 * (dz0 - dz1))};

    // If points are colinear, get the extrema and use the midpoint as the centre
    if (std::round(std::abs(scale)) == 0.f)
    {
        m_ccx = 0.5f * (minx + maxx);
        m_ccz = 0.5f * (minz + maxz);
        const float dx{m_ccx - minx};
        const float dz{m_ccz - minz};
        m_ccr2 = dx * dx + dz * dz;
    }
    else
    {
        m_ccx = ox + ((dx0 * dx0 + dz0 * dz0) * (dz1 - dz2) + (dx1 * dx1 + dz1 * dz1) * (dz2 - dz0) + (dx2 * dx2 + dz2 * dz2) * (dz0 - dz1)) / scale;
        m_ccz = oz + ((dx0 * dx0 + dz0 * dz0) * (dx2 - dx1) + (dx1 * dx1 + dz1 * dz1) * (dx0 - dx2) + (dx2 * dx2 + dz2 * dz2) * (dx1 - dx0)) / scale;
        const float dx{m_ccx - this->m_v0->m_x};
        const float dz{m_ccz - this->m_v0->m_z};
        m_ccr2 = dx * dx + dz * dz;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDelaunayTriangulationHelper::Triangle::GetEdges(EdgeVector &edges) const
{
    edges.emplace_back(this->m_e01);
    edges.emplace_back(this->m_e12);
    edges.emplace_back(this->m_e20);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDelaunayTriangulationHelper::Triangle::GetUniqueEdges(EdgeVector &edges) const
{
    if (std::find_if(edges.begin(), edges.end(), [this](const Edge *const &edge) { return *edge == *this->m_e01; }) == edges.end())
        edges.emplace_back(this->m_e01);
    if (std::find_if(edges.begin(), edges.end(), [this](const Edge *const &edge) { return *edge == *this->m_e12; }) == edges.end())
        edges.emplace_back(this->m_e12);
    if (std::find_if(edges.begin(), edges.end(), [this](const Edge *const &edge) { return *edge == *this->m_e20; }) == edges.end())
        edges.emplace_back(this->m_e20);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDelaunayTriangulationHelper::Triangle::SharesEdge(const Triangle &other) const
{
    return *(this->m_e01) == *(other.m_e01) || *(this->m_e01) == *(other.m_e12) || *(this->m_e01) == *(other.m_e20) ||
        *(this->m_e12) == *(other.m_e01) || *(this->m_e12) == *(other.m_e12) || *(this->m_e12) == *(other.m_e20) ||
        *(this->m_e20) == *(other.m_e01) || *(this->m_e20) == *(other.m_e12) || *(this->m_e20) == *(other.m_e20);
}

} // namespace lar_content
