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

        inline float DistanceSquared(const Vertex &other) const
        {
            return (this->m_x - other.m_x)*(this->m_x - other.m_x) + (this->m_z - other.m_z)*(this->m_z - other.m_z);
        };

        inline bool operator==(const Vertex &other) const { return this->m_x == other.m_x && this->m_z == other.m_z; };

        const float m_x;    ///< The x coordinate of the vertex
        const float m_z;    ///< The z coordinate of the vertex
        const pandora::CaloHit *const m_pCaloHit; ///< The calo hit represented by the vertex
    };
    typedef std::vector<const Vertex *> VertexVector;

    class Circle
    {
    public:
        /**
         *  @brief Constructor
         *
         *  @param  vertex0 An endpoint of the edge
         *  @param  vertex1 An endpoint of the edge
         **/
        Circle(const float x, const float z, const float r);

        /**
         *  @brief Checks if a vertex is contained within this circle
         *
         *  @param  vertex The vertex whose containment should be checked
         *
         *  @return true if the vertex is contained within the circle, false otherwise
         */
        bool Contains(const Vertex &vertex) const;

        const float m_x; ///< x coordinate of the circle
        const float m_z; ///< z coordinate of the circle
        const float m_r; ///< radius of the circle
        const float m_r2; ///< square of the radius of the circle
    };

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

        inline float LengthSquared() const
        {
            return (this->m_v0->m_x - this->m_v1->m_x)*(this->m_v0->m_x - this->m_v1->m_x) +
                (this->m_v0->m_z - this->m_v1->m_z)*(this->m_v0->m_z - this->m_v1->m_z);
        };


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
         *  @brief Checks if a vertex is within this triangle's circumcircle
         *
         *  @param  vertex The vertex to be checked for containment
         *
         *  @return Whether or not the vertex is contained within the circumcircle of this triangle
         */
        bool IsInCircumcircle(const Vertex &vertex) const;

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
        const Circle m_circle;    ///< The circumcircle associated with this triangle
    };
    typedef std::vector<const Triangle *> TriangleVector;

    /**
     *  @brief  Construct a circumcircle from a set of vertices
     *
     *  @param  v0 a vertex of the triangle from which the circle should be constructed
     *  @param  v1 a vertex of the triangle from which the circle should be constructed
     *  @param  v2 a vertex of the triangle from which the circle should be constructed
     *
     *  @return the circumcircle through the vertices
     */
    static Circle CalculateCircumcircle(const Vertex &v0, const Vertex &v1, const Vertex &v2);

    /**
     *  @brief  Determine the centre of a circumcircle
     *
     *  @param  dx0 the x distance between vertex 0 and vertex 1
     *  @param  dz0 the z distance between vertex 0 and vertex 1
     *  @param  dx1 the x distance between vertex 0 and vertex 2
     *  @param  dz1 the z distance between vertex 0 and vertex 2
     *
     *  @return the centre of the circumcircle
     */
    static Vertex GetCircumcircleCentre(const float dx0, const float dz0, const float dx1, const float dz1);

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

    /**
     *  @brief  Construct a Delaunay triangulation from a list of calo hits. Note the caller takes ownership of the vertices and triangles
     *          and is responsible for deleting them subsequently.
     *
     *  @param  caloHitList the calo hit list containing the hits to triangulate
     *  @param  vertices the vector of vertices to populate
     *  @param  triangles the vector of triangles to populate
     */
    static void Triangulate(const pandora::CaloHitList &caloHitList, VertexVector &vertices, TriangleVector &triangles);

private:
    /**
     *  @brief  Construct the minimum bounding circle for a set of vertices using Welzl's algorithm
     *
     *  @param  vertices the vertices to be bound
     *
     *  @return the circle bounding all vertices
     */
    static Circle Welzl(const VertexVector &vertices);

    static Circle WelzlRecursive(VertexVector &vertices, const int i, VertexVector boundary);

    /**
     *  @brief  Construct a circle from a set of boundary points
     *
     *  @param  boundary the points on the boundary of the circle to construct
     *
     *  @return the circle lying on the vertices
     */
    static Circle MakeCircle(const VertexVector &boundary);
};

} // namespace lar_content

#endif // #ifndef LAR_DELAUNAY_TRIANGULATION_HELPER_H
