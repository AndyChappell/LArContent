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
#include <iterator>
#include <limits>
#include <random>

using namespace pandora;

namespace lar_content
{

LArDelaunayTriangulationHelper::Triangle *LArDelaunayTriangulationHelper::MakeInitialBoundingTriangle(const VertexVector &vertices)
{
    const Circle enclosing{LArDelaunayTriangulationHelper::Welzl(vertices)};
    const float alpha{0};
    const float beta{2 * M_PI / 3};
    const float gamma{4 * M_PI / 3};
    const Vertex *const v0{new Vertex(enclosing.m_x + 2 * enclosing.m_r * std::cos(alpha), enclosing.m_z + 2 * enclosing.m_r * std::sin(alpha))};
    const Vertex *const v1{new Vertex(enclosing.m_x + 2 * enclosing.m_r * std::cos(beta), enclosing.m_z + 2 * enclosing.m_r * std::sin(beta))};
    const Vertex *const v2{new Vertex(enclosing.m_x + 2 * enclosing.m_r * std::cos(gamma), enclosing.m_z + 2 * enclosing.m_r * std::sin(gamma))};

    return new Triangle(v0, v1, v2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Circle LArDelaunayTriangulationHelper::Welzl(const VertexVector &vertices)
{
    VertexVector verticesCopy(vertices);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::shuffle(verticesCopy.begin(), verticesCopy.end(), rng);
    VertexVector boundary;
    return LArDelaunayTriangulationHelper::WelzlRecursive(verticesCopy, verticesCopy.size(), boundary);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Circle LArDelaunayTriangulationHelper::WelzlRecursive(VertexVector &vertices, const int n, VertexVector boundary)
{
    if (n == 0 || boundary.size() == 3)
        return LArDelaunayTriangulationHelper::MakeCircle(boundary);

    const int idx{rand() % n};
    const Vertex *const p{vertices[idx]};
    std::swap(vertices[idx], vertices[n - 1]);

    Circle circle(LArDelaunayTriangulationHelper::WelzlRecursive(vertices, n - 1, boundary));
    if (circle.Contains(*p))
        return circle;

    boundary.push_back(p);
    return LArDelaunayTriangulationHelper::WelzlRecursive(vertices, n - 1, boundary);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Circle LArDelaunayTriangulationHelper::MakeCircle(const VertexVector &boundary)
{
    const size_t n{boundary.size()};
    if (n == 0)
    {
        return Circle(0, 0, 0);
    }
    else if (n == 1)
    {
        const Vertex &v0{*boundary[0]};
        return Circle(v0.m_x, v0.m_z, 0);
    }
    else if (n == 2)
    {
        const Vertex &v0{*boundary[0]};
        const Vertex &v1{*boundary[1]};

        const float x{0.5f * (v0.m_x + v1.m_x)};
        const float z{0.5f * (v0.m_z + v1.m_z)};
        const float dx{x - v0.m_x};
        const float dz{z - v0.m_z};
        const float r{std::sqrt(dx * dx + dz * dz)};

        return Circle(x, z, r);
    }
    else
    {
        const Vertex &v0{*boundary[0]};
        const Vertex &v1{*boundary[1]};
        const Vertex &v2{*boundary[2]};

        return LArDelaunayTriangulationHelper::CalculateCircumcircle(v0, v1, v2);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDelaunayTriangulationHelper::AddVertex(const Vertex *const pVertex, TriangleVector &triangles)
{
    EdgeVector edges;
    // Remove triangles whose circumcircles contain the vertex and collect the unique edges
    for (auto iter = triangles.begin(); iter != triangles.end(); )
    {
        const Triangle *const triangle{*iter};
        if (triangle->IsInCircumcircle(*pVertex))
        {
            triangle->GetUniqueEdges(edges);
            iter = triangles.erase(iter);
            delete triangle;
        }
        else
        {
            ++iter;
        }
    }

    // Construct new triangles from the vertex and the unique edges
    for (const Edge *const pEdge : edges)
    {
        const Triangle *const triangle{new Triangle(pEdge->m_v0, pEdge->m_v1, pVertex)};

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

void LArDelaunayTriangulationHelper::Triangulate(const CaloHitList &caloHitList, VertexVector &vertices, TriangleVector &triangles)
{
    for (const CaloHit *pCaloHit : caloHitList)
    {
        vertices.emplace_back(new Vertex(pCaloHit));
    }

    const Triangle *bounds{MakeInitialBoundingTriangle(vertices)};
    triangles.emplace_back(bounds);
    for (const Vertex *const pVertex : vertices)
    {
        AddVertex(pVertex, triangles);
    }

    ShrinkWrap(*bounds, triangles);
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

LArDelaunayTriangulationHelper::Circle::Circle(const float x, const float z, const float r) :
    m_x{x},
    m_z{z},
    m_r{r},
    m_r2{r * r}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDelaunayTriangulationHelper::Circle::Contains(const Vertex &vertex) const
{
    const float dx{vertex.m_x - this->m_x};
    const float dz{vertex.m_z - this->m_z};

    return (dx * dx + dz * dz) <= this->m_r2;
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
    m_circle{LArDelaunayTriangulationHelper::CalculateCircumcircle(*this->m_v0, *this->m_v1, *this->m_v2)}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDelaunayTriangulationHelper::Triangle::IsInCircumcircle(const Vertex &vertex) const
{
    return this->m_circle.Contains(vertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Circle LArDelaunayTriangulationHelper::CalculateCircumcircle(const Vertex &v0, const Vertex &v1, const Vertex &v2)
{
    const LArDelaunayTriangulationHelper::Vertex &temp{LArDelaunayTriangulationHelper::GetCircumcircleCentre(v1.m_x - v0.m_x, v1.m_z - v0.m_z,
        v2.m_x - v0.m_x, v2.m_z - v0.m_z)};
    const LArDelaunayTriangulationHelper::Vertex centre(temp.m_x + v0.m_x, temp.m_z + v0.m_z);
    const float dx{centre.m_x - v0.m_x};
    const float dz{centre.m_z - v0.m_z};
    
    return LArDelaunayTriangulationHelper::Circle(centre.m_x, centre.m_z, std::sqrt(dx * dx + dz * dz));
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArDelaunayTriangulationHelper::Vertex LArDelaunayTriangulationHelper::GetCircumcircleCentre(const float dx0, const float dz0, const float dx1,
    const float dz1)
{
    const float d0{dx0 * dx0 + dz0 * dz0};
    const float d1{dx1 * dx1 + dz1 * dz1};
    const float d2{dx0 * dz1 - dz0 * dx1};

    return LArDelaunayTriangulationHelper::Vertex((dz1 * d0 - dz0 * d1) / (2 * d2), (dx0 * d1 - dx1 * d0) / (2 * d2));
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
