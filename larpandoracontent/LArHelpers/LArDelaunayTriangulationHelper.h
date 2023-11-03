/**
 *  @file   larpandoracontent/LArHelpers/LArDelaunayTriangulationHelper.h
 *
 *  @brief  Header file for the Delaunay triangulation helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DELAUNAY_TRIANGULATION_HELPER_H
#define LAR_DELAUNAY_TRIANGULATION_HELPER_H 1

#include "Objects/Cluster.h"

namespace lar_content
{

/**
 *  @brief  LArDelaunayTriangulationHelper class
 */
class LArDelaunayTriangulationHelper
{
public:
    class Vertex
    {
    public:
        /**
         *  @brief Constructor
         *
         *  @param  caloHit The calo hit to be represented by this vertex
         **/
        Vertex(const pandora::CaloHit *const pCaloHit);

        /**
         *  @brief Constructor for vertices of a bounding triangle (which has no associated CaloHits)
         *
         *  @param  x the x coordinate of the vertex
         *  @param  z the z coordinate of the vertex
         **/
        Vertex(const float x, const float z);

        inline bool operator==(const Vertex &other) const { return this->m_x == other.m_x && this->m_z == other.m_z; };

        const float m_x;    ///< The x coordinate of the vertex
        const float m_z;    ///< The z coordinate of the vertex
        const pandora::CaloHit *const m_pCaloHit; ///< The calo hit represented by the vertex
    };
    typedef std::vector<const Vertex *> VertexVector;

    class Edge
    {
    public:
        /**
         *  @brief Constructor
         *
         *  @param  vertex0 An endpoint of the edge
         *  @param  vertex1 An endpoint of the edge
         **/
        Edge(const Vertex *const pVertex0, const Vertex *const pVertex1);

        inline bool operator==(const Edge &other) const { return (*(this->m_v0) == *(other.m_v0) && *(this->m_v1) == *(other.m_v1)) ||
            (*(this->m_v0) == *(other.m_v1) && *(this->m_v1) == *(other.m_v0)); };

        const Vertex *const m_v0; ///< An endpoint of the edge
        const Vertex *const m_v1; ///< An endpoint of the edge
    };
    typedef std::vector<const Edge *> EdgeVector;

    class Triangle
    {
    public:
        /**
         *  @brief Constructor
         *
         *  @param  vertex0 A vertex of the triangle
         *  @param  vertex1 A vertex of the triangle
         *  @param  vertex2 A vertex of the triangle
         **/
        Triangle(const Vertex *const pVertex0, const Vertex *const pVertex1, const Vertex *const pVertex2);

        /**
         *  @brief Checks if a vertex is contained within this triangle's circumcircle
         *
         *  @param  vertex The vertex whose containment should be checked
         *
         *  @return true if the vertex is contained within the triangle's circumcircle, false otherwise
         */
        bool IsInCircumCircle(const Vertex *const vertex) const;

        /**
         *  @brief Retrieves the edges in this triangle
         *
         *  @param  edges The output vector of edges
         */
        void GetEdges(EdgeVector &edges) const;

        /**
         *  @brief Retrieves the edges in this triangle not already present in the vector
         *
         *  @param  edges The output set of edges
         */
        void GetUniqueEdges(EdgeVector &edges) const;

        /**
         *  @brief Checks if this triangle shares an edge with the specified triangle
         *
         *  @param  other The triangle to compare
         *
         *  @return true if this triangle shares an edge with the specified triangle, false otherwise
         */
        bool SharesEdge(const Triangle &other) const;

        const Vertex *const m_v0; ///< A vertex of the triangle
        const Vertex *const m_v1; ///< A vertex of the triangle
        const Vertex *const m_v2; ///< A vertex of the triangle
        const Edge *const m_e01; ///< The edge between vertex 0 and vertex 1 of the triangle
        const Edge *const m_e12; ///< The edge between vertex 1 and vertex 2 of the triangle
        const Edge *const m_e20; ///< The edge between vertex 2 and vertex 0 of the triangle

    private:
        void CalculateCircumCircle();

        float m_ccx;    ///< x coordinate of the circumcircle
        float m_ccz;    ///< z coordinate of the circumcircle
        float m_ccr2;   ///< square of the radius of the circumcircle
    };
    typedef std::vector<const Triangle *> TriangleVector;

    /**
     *  @brief  Construct the bounding triangle for a set of vertices
     *
     *  @param  vertices the vertices to be bound
     *
     *  @return the triangle bounding all vertices
     */
    static Triangle *MakeInitialBoundingTriangle(const VertexVector &vertices);

    /**
     *  @brief  Add a vertex to the set of triangles
     *
     *  @param  vertex the vertex to be added
     *  @param  triangles the list of triangles to be augmented by the new vertex
     */
    static void AddVertex(const Vertex *const pVertex, TriangleVector &triangles);

    /**
     *  @brief  Remove triangles that share an edge with the bounding triangle
     *
     *  @param  bounds the original bounding triangle
     *  @param  triangles the vector of triangles to be modified
     */
    static void ShrinkWrap(const Triangle &bounds, TriangleVector &triangles);
};

} // namespace lar_content

#endif // #ifndef LAR_DELAUNAY_TRIANGULATION_HELPER_H
