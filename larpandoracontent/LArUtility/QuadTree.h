/**
 *  @file   larpandoracontent/LArUtility/QuadTree.h
 *
 *  @brief  Header file for the quad tree class
 *
 *  $Log: $
 */
#ifndef LAR_QUAD_TREE_H
#define LAR_QUAD_TREE_H 1

#include "Objects/CaloHit.h"
#include "Objects/CartesianVector.h"

#include "Pandora/PandoraInternal.h"

#include <array>
#include <vector>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  QuadTree defines a structured bounding region.
 */
class QuadTree
{
public:
    enum directions{LEFT=0, RIGHT, UPSTREAM=0, DOWNSTREAM};
    /**
     *  @brief  Construct a QuadTree from a given list of calo hits
     *
     *  @param  caloHitList The list of calo hits to be structured
     */
    QuadTree(const pandora::CaloHitList &caloHitList);

    ~QuadTree();

    /**
     *  @brief  Collect all of the hits within the given bounds
     *
     *  @param  xMin The minimum x-coordinate for collecting a hit
     *  @param  zMin The minimum z-coordinate for collecting a hit
     *  @param  xMax The maximum x-coordinate for collecting a hit
     *  @param  zMax The maximum z-coordinate for collecting a hit
     *  @param  backgroundHitList The output list of collected hits
     */
    void Find(const float xMin, const float zMin, const float xMax, const float zMax, pandora::CaloHitList &backgroundHitList) const;

    /**
     *  @brief  Produces a string representation of the tree structure
     *
     *  @return A string representation of the tree
     */
    const std::string ToString() const;

    /**
     *  @brief  A node describing a sub-region within the quad tree
     */
    class Node
    {
    public:
        /**
         *  @brief  Construct a Node from a given list of calo hits with a specified set of bounds
         *
         *  @param  caloHitList The list of calo hits to be structured
         *  @param  xMin The inclusive lower bound for the x-coordinate
         *  @param  xMax The exclusive upper bound for the x-coordinate
         *  @param  zMin The inclusive lower bound for the x-coordinate
         *  @param  zMax The exclusive upper bound for the x-coordinate
         */
        Node(const pandora::CaloHitList &caloHitList, const float xMin, const float xMax, const float zMin, const float zMax);

        ~Node();

        /**
         *  @brief  Collect all of the hits within the given bounds
         *
         *  @param  xMin The minimum x-coordinate for collecting a hit
         *  @param  zMin The minimum z-coordinate for collecting a hit
         *  @param  xMax The maximum x-coordinate for collecting a hit
         *  @param  zMax The maximum z-coordinate for collecting a hit
         *  @param  backgroundHitList The output list of collected hits
         */
        void Find(const float xMin, const float zMin, const float xMax, const float zMax, pandora::CaloHitList &backgroundHitList) const;

        /**
         *  Retrieve the position of this node
         *
         *  @return  The position of this node
         */
        const pandora::CartesianVector &GetPosition() const;

        /**
         *  Retrieve the list of hits in this node
         *
         *  @return The list of hits in this node
         */
        const pandora::CaloHitList &GetCaloHitList() const;

        /**
         *  Retrieve the child node representing the specified region
         *
         *  @param  xRegion The x region of interest (LEFT, RIGHT)
         *  @param  zRegion The z region of interest (UPSTREAM, DOWNSTREAM)
         *  @return The node representing the region of interest, or nullptr if this region has no hits
         */
        const Node *GetNode(const enum directions xRegion, const enum directions zRegion) const;

        /**
         *  @brief  Produces a string representation of the tree structure
         *
         *  @return A string representation of the tree
         */
        const std::string ToString(const std::string &tab = "") const;

    protected:
        Node();

    private:
        pandora::CartesianVector m_centre;  ///< The centre of the bounding region
        pandora::CartesianVector m_min;     ///< The minimum coordinate of the bounding region
        pandora::CartesianVector m_max;     ///< The maximum coordinate of the bounding region
        std::array<std::array<const Node *, 2>, 2> m_children;
        pandora::CaloHitList m_caloHitList;          ///< The calo hits contained within the bounding region
    };

private:
    Node *m_root;
};

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // LAR_QUAD_TREE_H
