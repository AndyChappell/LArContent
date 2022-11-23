/**
 *  @file   larpandoracontent/LArUtility/QuadTree.cc
 *
 *  @brief  Implementation of the quad tree class
 *
 *  $Log: $
 */

#include "larpandoracontent/LArUtility/QuadTree.h"

using namespace pandora;

namespace lar_content
{

QuadTree::QuadTree(const CaloHitList &caloHitList) :
    m_root{nullptr}
{
    float xMin{std::numeric_limits<float>::max()}, xMax{-std::numeric_limits<float>::max()};
    float zMin{std::numeric_limits<float>::max()}, zMax{-std::numeric_limits<float>::max()};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        xMin = std::min(xMin, pos.GetX() - 1);
        xMax = std::max(xMax, pos.GetX() + 1);
        zMin = std::min(zMin, pos.GetZ() - 1);
        zMax = std::max(zMax, pos.GetZ() + 1);
    }

    m_root = new Node(caloHitList, xMin, xMax, zMin, zMax);
}

//------------------------------------------------------------------------------------------------------------------------------------------

QuadTree::~QuadTree()
{
    if (m_root)
        delete m_root;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void QuadTree::Find(const float xMin, const float zMin, const float xMax, const float zMax, CaloHitList &backgroundHitList) const
{
    m_root->Find(xMin, zMin, xMax, zMax, backgroundHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string QuadTree::ToString() const
{
    return m_root->ToString();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

QuadTree::Node::Node(const CaloHitList &caloHitList, const float xMin, const float xMax, const float zMin, const float zMax) :
    m_centre{CartesianVector(0.5f * (xMin + xMax), 0, 0.5f * (zMin + zMax))},
    m_min{CartesianVector(xMin, 0, zMin)},
    m_max{CartesianVector(xMax, 0, zMax)}
{
    for (enum directions x : {LEFT, RIGHT})
        for (enum directions z : {UPSTREAM, DOWNSTREAM})
            m_children[x][z] = nullptr;

    if (((xMax - xMin) > 5.f) && ((zMax - zMin) > 5.f))
    {
        CaloHitList hitListLU, hitListLD, hitListRU, hitListRD;
        for (const CaloHit *pCaloHit : caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            if (m_min.GetX() <= pos.GetX() && pos.GetX() < m_centre.GetX())
            {
                if (m_min.GetZ() <= pos.GetZ() && pos.GetZ() < m_centre.GetZ())
                    hitListLU.emplace_back(pCaloHit);
                else if (m_centre.GetZ() <= pos.GetZ() && pos.GetZ() < m_max.GetZ())
                    hitListLD.emplace_back(pCaloHit);
            }
            else if (m_centre.GetX() <= pos.GetX() && pos.GetX() < m_max.GetX())
            {
                if (m_min.GetZ() <= pos.GetZ() && pos.GetZ() < m_centre.GetZ())
                    hitListRU.emplace_back(pCaloHit);
                else if (m_centre.GetZ() <= pos.GetZ() && pos.GetZ() < m_max.GetZ())
                    hitListRD.emplace_back(pCaloHit);
            }
        }

        if (!hitListLU.empty())
            m_children[LEFT][UPSTREAM] = new Node(hitListLU, m_min.GetX(), m_centre.GetX(), m_min.GetZ(), m_centre.GetZ());
        if (!hitListLD.empty())
            m_children[LEFT][DOWNSTREAM] = new Node(hitListLD, m_min.GetX(), m_centre.GetX(), m_centre.GetZ(), m_max.GetZ());
        if (!hitListRU.empty())
            m_children[RIGHT][UPSTREAM] = new Node(hitListRU, m_centre.GetX(), m_max.GetX(), m_min.GetZ(), m_centre.GetZ());
        if (!hitListRD.empty())
            m_children[RIGHT][DOWNSTREAM] = new Node(hitListRD, m_centre.GetX(), m_max.GetX(), m_centre.GetZ(), m_max.GetZ());
    }
    else
    {
        for (const CaloHit *const pCaloHit : caloHitList)
            m_caloHitList.emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

QuadTree::Node::~Node()
{
    for (auto x : {LEFT, RIGHT})
        for (auto z : {UPSTREAM, DOWNSTREAM})
            if (m_children[x][z])
                delete m_children[x][z];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void QuadTree::Node::Find(const float xMin, const float zMin, const float xMax, const float zMax, CaloHitList &backgroundHitList) const
{
    if (!m_caloHitList.empty())
    {
        for (const CaloHit *const pCaloHit : m_caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            if (xMin <= pos.GetX() && pos.GetX() <= xMax && zMin <= pos.GetZ() && pos.GetZ() <= zMax)
                backgroundHitList.emplace_back(pCaloHit);
        }
    }
    else
    {
        if (xMin <= m_centre.GetX())
        {
            if (zMin <= m_centre.GetZ())
            {
                if (m_children[LEFT][UPSTREAM])
                    m_children[LEFT][UPSTREAM]->Find(xMin, zMin, xMax, zMax, backgroundHitList);
            }
            if (m_centre.GetZ() <= zMax)
            {
                if (m_children[LEFT][DOWNSTREAM])
                    m_children[LEFT][DOWNSTREAM]->Find(xMin, zMin, xMax, zMax, backgroundHitList);
            }
        }
        if (m_centre.GetX() <= xMax)
        {
            if (zMin <= m_centre.GetZ())
            {
                if (m_children[RIGHT][UPSTREAM])
                    m_children[RIGHT][UPSTREAM]->Find(xMin, zMin, xMax, zMax, backgroundHitList);
            }
            if (m_centre.GetZ() <= zMax)
            {
                if (m_children[RIGHT][DOWNSTREAM])
                    m_children[RIGHT][DOWNSTREAM]->Find(xMin, zMin, xMax, zMax, backgroundHitList);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &QuadTree::Node::GetPosition() const
{
    return m_centre;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &QuadTree::Node::GetCaloHitList() const
{
    return m_caloHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const QuadTree::Node *QuadTree::Node::GetNode(const enum directions xRegion, const enum directions zRegion) const
{
    if ((xRegion != LEFT && xRegion != RIGHT) || (zRegion != UPSTREAM && zRegion != DOWNSTREAM))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return m_children[xRegion][zRegion];
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string QuadTree::Node::ToString(const std::string &tab) const
{
    std::string str{tab + "Centre [" + std::to_string(m_centre.GetX()) + "," + std::to_string(m_centre.GetZ()) + "]\n"};
    for (enum directions x : { LEFT, RIGHT })
    {
        for (enum directions z : { UPSTREAM, DOWNSTREAM })
        {
            if (m_children[x][z])
            {
                str += m_children[x][z]->ToString(tab + "  ") + "\n";
            }
        }
    }
    if (!m_caloHitList.empty())
        str += tab + std::string("  ") + std::to_string(m_caloHitList.size()) + "hits\n";

    return str;
}

} // namespace lar_content

