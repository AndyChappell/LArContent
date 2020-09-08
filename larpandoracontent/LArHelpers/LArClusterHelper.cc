/**
 *  @file   larpandoracontent/LArHelpers/LArClusterHelper.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <list>

using namespace pandora;

namespace lar_content
{

HitType LArClusterHelper::GetClusterHitType(const Cluster *const pCluster)
{
    if (0 == pCluster->GetNCaloHits())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // TODO Find a way of confirming single hit-type clustering mode; currently only checked in ListPreparation algorithm
    // if (!pandora->GetSettings()->SingleHitTypeClusteringMode())
    //     throw StatusCodeException(STATUS_CODE_FAILURE);

    return (*(pCluster->GetOrderedCaloHitList().begin()->second->begin()))->GetHitType();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClustersUVW(const ClusterList &inputClusters, ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW)
{
    for (ClusterList::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(*iter));

        if (TPC_VIEW_U == hitType) clusterListU.push_back(*iter);
        else if (TPC_VIEW_V == hitType) clusterListV.push_back(*iter);
        else if (TPC_VIEW_W == hitType) clusterListW.push_back(*iter);
        else throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClustersByHitType(const ClusterList &inputClusters, const HitType hitType, ClusterList &clusterList)
{
    for (ClusterList::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        if (hitType == LArClusterHelper::GetClusterHitType(*iter))
            clusterList.push_back(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLengthSquared(const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // ATTN In 2D case, we will actually calculate the quadrature sum of deltaX and deltaU/V/W
    float minX(std::numeric_limits<float>::max()), maxX(-std::numeric_limits<float>::max());
    float minY(std::numeric_limits<float>::max()), maxY(-std::numeric_limits<float>::max());
    float minZ(std::numeric_limits<float>::max()), maxZ(-std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CartesianVector &hitPosition((*hitIter)->GetPositionVector());
            minX = std::min(hitPosition.GetX(), minX);
            maxX = std::max(hitPosition.GetX(), maxX);
            minY = std::min(hitPosition.GetY(), minY);
            maxY = std::max(hitPosition.GetY(), maxY);
            minZ = std::min(hitPosition.GetZ(), minZ);
            maxZ = std::max(hitPosition.GetZ(), maxZ);
        }
    }

    const float deltaX(maxX - minX), deltaY(maxY - minY), deltaZ(maxZ - minZ);
    return (deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLength(const Cluster *const pCluster)
{
    return std::sqrt(LArClusterHelper::GetLengthSquared(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetEnergyFromLength(const Cluster *const pCluster)
{
    const float dEdX(0.002f); // approximately 2 MeV/cm
    return (dEdX * LArClusterHelper::GetLength(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArClusterHelper::GetLayerSpan(const Cluster *const pCluster)
{
    return (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster)
{
    const unsigned int nOccupiedLayers(pCluster->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    const unsigned int nOccupiedLayers(pCluster1->GetOrderedCaloHitList().size() + pCluster2->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + std::max(pCluster1->GetOuterPseudoLayer(), pCluster2->GetOuterPseudoLayer()) -
        std::min(pCluster1->GetInnerPseudoLayer(), pCluster2->GetInnerPseudoLayer()));

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const ClusterList &clusterList1, const ClusterList &clusterList2)
{
    if (clusterList1.empty() || clusterList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float closestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iterEnd1 = clusterList1.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *const pCluster1 = *iter1;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster1, clusterList2));

        if (thisDistance < closestDistance)
            closestDistance = thisDistance;
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const Cluster *const pCluster, const ClusterList &clusterList)
{
    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float closestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pTestCluster = *iter;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster, pTestCluster));

        if (thisDistance < closestDistance)
            closestDistance = thisDistance;
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    CartesianVector closestPosition1(0.f, 0.f, 0.f);
    CartesianVector closestPosition2(0.f, 0.f, 0.f);

    LArClusterHelper::GetClosestPositions(pCluster1, pCluster2, closestPosition1, closestPosition2);

    return (closestPosition1 - closestPosition2).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const ClusterList &clusterList)
{
    return (position - LArClusterHelper::GetClosestPosition(position, clusterList)).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const Cluster *const pCluster)
{
    return (position - LArClusterHelper::GetClosestPosition(position, pCluster)).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::GetClosestPosition(const CartesianVector &position, const ClusterList &clusterList)
{
    bool distanceFound(false);
    float closestDistanceSquared(std::numeric_limits<float>::max());
    CartesianVector closestPosition(0.f, 0.f, 0.f);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pTestCluster = *iter;
        const CartesianVector thisPosition(LArClusterHelper::GetClosestPosition(position, pTestCluster));
        const float thisDistanceSquared((position - thisPosition).GetMagnitudeSquared());

        if (thisDistanceSquared < closestDistanceSquared)
        {
            distanceFound = true;
            closestDistanceSquared = thisDistanceSquared;
            closestPosition = thisPosition;
        }
    }

    if (distanceFound)
        return closestPosition;

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::GetClosestPosition(const CartesianVector &position, const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    const CaloHit *pClosestCaloHit(NULL);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;
            const float distanceSquared((pCaloHit->GetPositionVector() - position).GetMagnitudeSquared());

            if (distanceSquared < closestDistanceSquared)
            {
                closestDistanceSquared = distanceSquared;
                pClosestCaloHit = pCaloHit;
            }
        }
    }

    if (pClosestCaloHit)
        return pClosestCaloHit->GetPositionVector();

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClosestPositions(const Cluster *const pCluster1, const Cluster *const pCluster2, CartesianVector &outputPosition1,
    CartesianVector &outputPosition2)
{
    bool distanceFound(false);
    float minDistanceSquared(std::numeric_limits<float>::max());

    CartesianVector closestPosition1(0.f, 0.f, 0.f);
    CartesianVector closestPosition2(0.f, 0.f, 0.f);

    const OrderedCaloHitList &orderedCaloHitList1(pCluster1->GetOrderedCaloHitList());
    const OrderedCaloHitList &orderedCaloHitList2(pCluster2->GetOrderedCaloHitList());

    // Loop over hits in cluster 1
    for (OrderedCaloHitList::const_iterator iter1 = orderedCaloHitList1.begin(), iter1End = orderedCaloHitList1.end(); iter1 != iter1End; ++iter1)
    {
        for (CaloHitList::const_iterator hitIter1 = iter1->second->begin(), hitIter1End = iter1->second->end(); hitIter1 != hitIter1End; ++hitIter1)
        {
            const CartesianVector &positionVector1((*hitIter1)->GetPositionVector());

            // Loop over hits in cluster 2
            for (OrderedCaloHitList::const_iterator iter2 = orderedCaloHitList2.begin(), iter2End = orderedCaloHitList2.end(); iter2 != iter2End; ++iter2)
            {
                for (CaloHitList::const_iterator hitIter2 = iter2->second->begin(), hitIter2End = iter2->second->end(); hitIter2 != hitIter2End; ++hitIter2)
                {
                    const CartesianVector &positionVector2((*hitIter2)->GetPositionVector());

                    const float distanceSquared((positionVector1 - positionVector2).GetMagnitudeSquared());

                    if (distanceSquared < minDistanceSquared)
                    {
                        minDistanceSquared = distanceSquared;
                        closestPosition1 = positionVector1;
                        closestPosition2 = positionVector2;
                        distanceFound = true;
                    }
                }
            }
        }
    }

    if (!distanceFound)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    outputPosition1 = closestPosition1;
    outputPosition2 = closestPosition2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetClusterBoundingBox(const Cluster *const pCluster, CartesianVector &minimumCoordinate, CartesianVector &maximumCoordinate)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float xmin(std::numeric_limits<float>::max());
    float ymin(std::numeric_limits<float>::max());
    float zmin(std::numeric_limits<float>::max());
    float xmax(-std::numeric_limits<float>::max());
    float ymax(-std::numeric_limits<float>::max());
    float zmax(-std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());
            xmin = std::min(hit.GetX(), xmin);
            xmax = std::max(hit.GetX(), xmax);
            ymin = std::min(hit.GetY(), ymin);
            ymax = std::max(hit.GetY(), ymax);
            zmin = std::min(hit.GetZ(), zmin);
            zmax = std::max(hit.GetZ(), zmax);
        }
    }

    minimumCoordinate.SetValues(xmin, ymin, zmin);
    maximumCoordinate.SetValues(xmax, ymax, zmax);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArClusterHelper::GetAverageZ(const Cluster *const pCluster, const float xmin, const float xmax, float &averageZ)
{
    averageZ = std::numeric_limits<float>::max();

    if (xmin > xmax)
        return STATUS_CODE_INVALID_PARAMETER;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float zsum(0.f);
    int count(0);

    for (OrderedCaloHitList::const_iterator ochIter = orderedCaloHitList.begin(), ochIterEnd = orderedCaloHitList.end(); ochIter != ochIterEnd; ++ochIter)
    {
        for (CaloHitList::const_iterator hIter = ochIter->second->begin(), hIterEnd = ochIter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;
            const CartesianVector &hit(pCaloHit->GetPositionVector());

            if (hit.GetX() < xmin || hit.GetX() > xmax)
                continue;

            zsum += hit.GetZ();
            ++count;
        }
    }

    if (count == 0)
        return STATUS_CODE_NOT_FOUND;

    averageZ = zsum / static_cast<float>(count);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const ClusterList &clusterList, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate)
{
    OrderedCaloHitList orderedCaloHitList;

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(pCluster->GetOrderedCaloHitList()));
    }

    return LArClusterHelper::GetExtremalCoordinates(orderedCaloHitList, innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const Cluster *const pCluster, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate)
{
  return LArClusterHelper::GetExtremalCoordinates(pCluster->GetOrderedCaloHitList(), innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const OrderedCaloHitList &orderedCaloHitList, CartesianVector &innerCoordinate,
    CartesianVector &outerCoordinate)
{
    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianPointVector coordinateVector;

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;
            coordinateVector.push_back(pCaloHit->GetPositionVector());
        }
    }

    std::sort(coordinateVector.begin(), coordinateVector.end(), LArClusterHelper::SortCoordinatesByPosition);
    return LArClusterHelper::GetExtremalCoordinates(coordinateVector, innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetExtremalCoordinates(const CartesianPointVector &coordinateVector, CartesianVector &innerCoordinate,
    CartesianVector &outerCoordinate)
{
    if (coordinateVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Find the extremal values of the X, Y and Z coordinates
    float xMin(+std::numeric_limits<float>::max());
    float yMin(+std::numeric_limits<float>::max());
    float zMin(+std::numeric_limits<float>::max());
    float xMax(-std::numeric_limits<float>::max());
    float yMax(-std::numeric_limits<float>::max());
    float zMax(-std::numeric_limits<float>::max());

    for (CartesianPointVector::const_iterator pIter = coordinateVector.begin(), pIterEnd = coordinateVector.end(); pIter != pIterEnd; ++pIter)
    {
        const CartesianVector &pos = *pIter;
        xMin = std::min(pos.GetX(), xMin);
        xMax = std::max(pos.GetX(), xMax);
        yMin = std::min(pos.GetY(), yMin);
        yMax = std::max(pos.GetY(), yMax);
        zMin = std::min(pos.GetZ(), zMin);
        zMax = std::max(pos.GetZ(), zMax);
    }

    // Choose the coordinate with the greatest span (keeping any ties)
    const float xAve(0.5f * (xMin + xMax));
    const float yAve(0.5f * (yMin + yMax));
    const float zAve(0.5f * (zMin + zMax));

    const float xSpan(xMax - xMin);
    const float ySpan(yMax - yMin);
    const float zSpan(zMax - zMin);

    const bool useX((xSpan > std::numeric_limits<float>::epsilon()) && (xSpan + std::numeric_limits<float>::epsilon() > std::max(ySpan, zSpan)));
    const bool useY((ySpan > std::numeric_limits<float>::epsilon()) && (ySpan + std::numeric_limits<float>::epsilon() > std::max(zSpan, xSpan)));
    const bool useZ((zSpan > std::numeric_limits<float>::epsilon()) && (zSpan + std::numeric_limits<float>::epsilon() > std::max(xSpan, ySpan)));

    // Find the extremal hits separately for the chosen coordinates
    CartesianPointVector candidateVector;

    for (CartesianPointVector::const_iterator pIter = coordinateVector.begin(), pIterEnd = coordinateVector.end(); pIter != pIterEnd; ++pIter)
    {
        const CartesianVector &pos = *pIter;

        if (useX)
        {
            if (((pos.GetX() - xMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetX() - xMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }

        if (useY)
        {
            if (((pos.GetY() - yMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetY() - yMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }

        if (useZ)
        {
            if (((pos.GetZ() - zMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetZ() - zMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }
    }

    // Finally, find the pair of hits that are separated by the greatest distance
    CartesianVector firstCoordinate(xAve, yAve, zAve);
    CartesianVector secondCoordinate(xAve, yAve, zAve);
    float maxDistanceSquared(+std::numeric_limits<float>::epsilon());

    for (CartesianPointVector::const_iterator iterI = candidateVector.begin(), iterEndI = candidateVector.end(); iterI != iterEndI; ++iterI)
    {
        const CartesianVector &posI = *iterI;

        for (CartesianPointVector::const_iterator iterJ = iterI, iterEndJ = candidateVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            const CartesianVector &posJ = *iterJ;

            const float distanceSquared((posI - posJ).GetMagnitudeSquared());

            if (distanceSquared > maxDistanceSquared)
            {
                maxDistanceSquared = distanceSquared;
                firstCoordinate = posI;
                secondCoordinate = posJ;
            }
        }
    }

    // Set the inner and outer coordinates (Check Z first, then X in the event of a tie)
    const float deltaZ(secondCoordinate.GetZ() - firstCoordinate.GetZ());
    const float deltaX(secondCoordinate.GetX() - firstCoordinate.GetX());

    if ((deltaZ > 0.f) || ((std::fabs(deltaZ) < std::numeric_limits<float>::epsilon()) && (deltaX > 0.f)))
    {
        innerCoordinate = firstCoordinate;
        outerCoordinate = secondCoordinate;
    }
    else
    {
        innerCoordinate = secondCoordinate;
        outerCoordinate = firstCoordinate;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetCoordinateVector(const Cluster *const pCluster, CartesianPointVector &coordinateVector)
{
    for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
    {
        for (const CaloHit *const pCaloHit : *layerEntry.second)
            coordinateVector.push_back(pCaloHit->GetPositionVector());
    }

    std::sort(coordinateVector.begin(), coordinateVector.end(), LArClusterHelper::SortCoordinatesByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetCaloHitListInBoundingBox(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBound,
        const pandora::CartesianVector &upperBound, pandora::CaloHitList &caloHitList)
{
    const bool useX(std::fabs(upperBound.GetX() - lowerBound.GetX()) > std::numeric_limits<float>::epsilon());
    const bool useY(std::fabs(upperBound.GetY() - lowerBound.GetY()) > std::numeric_limits<float>::epsilon());
    const bool useZ(std::fabs(upperBound.GetZ() - lowerBound.GetZ()) > std::numeric_limits<float>::epsilon());
    if (!useX && !useY && !useZ)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float minX(std::min(lowerBound.GetX(), upperBound.GetX()));
    const float maxX(std::max(lowerBound.GetX(), upperBound.GetX()));
    const float minY(std::min(lowerBound.GetY(), upperBound.GetY()));
    const float maxY(std::max(lowerBound.GetY(), upperBound.GetY()));
    const float minZ(std::min(lowerBound.GetZ(), upperBound.GetZ()));
    const float maxZ(std::max(lowerBound.GetZ(), upperBound.GetZ()));

    for (const OrderedCaloHitList::value_type &layerEntry : pCluster->GetOrderedCaloHitList())
    {
        for (const CaloHit *const pCaloHit : *layerEntry.second)
        {
            const CartesianVector &hitPosition = pCaloHit->GetPositionVector();
            if (useX &&
                (hitPosition.GetX() < minX - std::numeric_limits<float>::epsilon() ||
                 hitPosition.GetX() > maxX + std::numeric_limits<float>::epsilon()))
                    continue;
            else if (useY &&
                (hitPosition.GetY() < minY - std::numeric_limits<float>::epsilon() ||
                 hitPosition.GetY() > maxY + std::numeric_limits<float>::epsilon()))
                    continue;
            else if (useZ &&
                (hitPosition.GetZ() < minZ - std::numeric_limits<float>::epsilon() ||
                 hitPosition.GetZ() > maxZ + std::numeric_limits<float>::epsilon()))
                    continue;

            caloHitList.push_back(pCaloHit);
        }
    }
    caloHitList.sort(LArClusterHelper::SortHitsByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNOccupiedLayers(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nOccupiedLayersLhs(pLhs->GetOrderedCaloHitList().size());
    const unsigned int nOccupiedLayersRhs(pRhs->GetOrderedCaloHitList().size());

    if (nOccupiedLayersLhs != nOccupiedLayersRhs)
        return (nOccupiedLayersLhs > nOccupiedLayersRhs);

    return SortByNHits(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNHits(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nHitsLhs(pLhs->GetNCaloHits());
    const unsigned int nHitsRhs(pRhs->GetNCaloHits());

    if (nHitsLhs != nHitsRhs)
        return (nHitsLhs > nHitsRhs);

    return SortByLayerSpan(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByLayerSpan(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return SortByInnerLayer(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByInnerLayer(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int innerLayerLhs(pLhs->GetInnerPseudoLayer());
    const unsigned int innerLayerRhs(pRhs->GetInnerPseudoLayer());

    if (innerLayerLhs != innerLayerRhs)
      return (innerLayerLhs < innerLayerRhs);

    return SortByPosition(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByPosition(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const CartesianVector deltaPositionIL(pRhs->GetCentroid(pRhs->GetInnerPseudoLayer()) - pLhs->GetCentroid(pLhs->GetInnerPseudoLayer()));

    if (std::fabs(deltaPositionIL.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionIL.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionIL.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionIL.GetX() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionIL.GetY()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionIL.GetY() > std::numeric_limits<float>::epsilon());

    const CartesianVector deltaPositionOL(pRhs->GetCentroid(pRhs->GetOuterPseudoLayer()) - pLhs->GetCentroid(pLhs->GetOuterPseudoLayer()));

    if (std::fabs(deltaPositionOL.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionOL.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionOL.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionOL.GetX() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPositionOL.GetY()) > std::numeric_limits<float>::epsilon())
        return (deltaPositionOL.GetY() > std::numeric_limits<float>::epsilon());

    // Use pulse height to resolve ties
    return SortByPulseHeight(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByPulseHeight(const Cluster *const pLhs, const Cluster *const pRhs)
{
    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortHitsByPosition(const CaloHit *const pLhs, const CaloHit *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPositionVector() - pLhs->GetPositionVector());

    if (std::fabs(deltaPosition.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetY()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetY() > std::numeric_limits<float>::epsilon());

    // Use pulse height to resolve ties
    return SortHitsByPulseHeight(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortHitsByPositionInX(const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPositionVector() - pLhs->GetPositionVector());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    return SortHitsByPosition(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortHitsByPulseHeight(const CaloHit *const pLhs, const CaloHit *const pRhs)
{
    // TODO: Think about the correct energy to use here (should they ever be different)
    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortCoordinatesByPosition(const CartesianVector &lhs, const CartesianVector &rhs)
{
    const CartesianVector deltaPosition(rhs - lhs);

    if (std::fabs(deltaPosition.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    return (deltaPosition.GetY() > std::numeric_limits<float>::epsilon());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetConcaveHull(const Cluster *const pCluster, const int k, CaloHitList &concaveHull)
{
    const OrderedCaloHitList &caloHits{pCluster->GetOrderedCaloHitList()};
    int nPoints{0};
    for (auto iter = caloHits.begin(); iter != caloHits.end(); ++iter)
        nPoints += iter->second->size();
    if (nPoints <= 10)
    {   // Strictly < 3 points doesn't form a concave hull, but we still want the points anyway
        // 10 hits or fewer is rarely wortwhile to compute, just keep everything
        for (auto iter = caloHits.begin(); iter != caloHits.end(); ++iter)
            for (auto hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter)
                concaveHull.push_back(*hitIter);
        return;
    }
    // In general k should be at least 3, but if there aren't at least k neighbours left in the dataset, set it to the number of neighbours
    const int kk{k >= 3 ?
        std::min(nPoints - 1, k) :
        std::min(nPoints - 1, 3)};
    // Set up the input hits - want them sorted, but then want a deque for efficient insertion/removal at front
    CaloHitList inputCaloHitList;
    for (auto iter = caloHits.begin(); iter != caloHits.end(); ++iter)
        for (auto hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter)
            inputCaloHitList.push_back(*hitIter);
    inputCaloHitList.sort(LArClusterHelper::SortHitsByPosition);
    std::deque<const CaloHit*> inputCaloHitQueue;
    for (const CaloHit *pCaloHit : inputCaloHitList)
        inputCaloHitQueue.push_back(pCaloHit);

    // Initialise the hull with the first hit (lowest z) and remove from the input queue
    const CaloHit *pFirstHit{inputCaloHitQueue.front()};
    const CaloHit *pCurrentHit{pFirstHit};
    inputCaloHitQueue.pop_front();
    concaveHull.push_back(pFirstHit);
    float xmin{pFirstHit->GetPositionVector().GetX()}, zmin{pFirstHit->GetPositionVector().GetZ()};
    float xmax{xmin}, zmax{zmin};

    int step{2};
    CartesianVector baseline{-1.f, 0.f, 0.f};
    typedef std::map<float, const CaloHit*, std::greater<float>> DistanceNeighbourMap;
    DistanceNeighbourMap distanceNeighbourMap;
    while ((pCurrentHit != pFirstHit || step == 2) && !inputCaloHitQueue.empty())
    {
        if (step == 5)
        {   // Need to add the starting point back to the candidate list to be able to close the hull, but must wait until 4 points have
            // been added or you trivially end up with a triangle
            inputCaloHitQueue.push_back(pFirstHit);
        }
        // Find the kk nearest neighbours to the current hit
        distanceNeighbourMap.clear();
        for (const CaloHit *pCaloHit : inputCaloHitQueue)
        {
            if (pCurrentHit == pCaloHit)
                continue;
            const float delta{(pCurrentHit->GetPositionVector() - pCaloHit->GetPositionVector()).GetMagnitudeSquared()};
            const int size{static_cast<int>(distanceNeighbourMap.size())};
            if (size < kk)
            {   // Don't yet have kk nearest neighbours, so just add it to the map
                distanceNeighbourMap.emplace(delta, pCaloHit);
            }
            else
            {   // Already have kk neighbours, check if the current hit is closer
                // Map is automatically sorted by descending distance, so only need to check the first element
                auto iter = distanceNeighbourMap.begin();
                if (delta < iter->first)
                {
                    distanceNeighbourMap.erase(iter);
                    distanceNeighbourMap.emplace(delta, pCaloHit);
                }
            }
        }
        if (step < 5 && distanceNeighbourMap.size() < kk)
        {   // Can get some edge cases where the first hit hasn't been added back, but there are too few neighbours
            const float delta{(pCurrentHit->GetPositionVector() - pFirstHit->GetPositionVector()).GetMagnitudeSquared()};
            distanceNeighbourMap.emplace(delta, pFirstHit);
        }

        // Calculated the clockwise angle between most recent edge and each candidate edge (default to horizontal baseline if current
        // hull has less than two hits
        DistanceNeighbourMap angleNeighbourMap;
        for (auto iter = distanceNeighbourMap.begin(); iter != distanceNeighbourMap.end(); ++iter)
        {
            // Order of vertices in baseline and new edge matters for a clockwise result
            // +z
            // ^
            // |
            //  ---> +x
            const CartesianVector &proposedVec{iter->second->GetPositionVector() - pCurrentHit->GetPositionVector()};
            float dot2d{proposedVec.GetX() * baseline.GetX() + proposedVec.GetZ() * baseline.GetZ()};
            float det2d{proposedVec.GetX() * baseline.GetZ() - proposedVec.GetZ() * baseline.GetX()};
            // Clockwise angle can be computed from atan2 by using fact that dot product is proportional to cosine, while determinant
            // is proportional to sine. Range is [-180, 180], so add 360 if < 0.
            float angle{std::atan2(det2d, dot2d)};
            if (angle < 0.f)
                angle += static_cast<float>(2 * M_PI);
            angleNeighbourMap.emplace(angle, iter->second);
        }
        if (angleNeighbourMap.empty())
            break;
        bool intersects{true};
        const CaloHit *pCandidate{nullptr};
        for (auto iter = angleNeighbourMap.begin(); iter != angleNeighbourMap.end() && intersects; ++iter)
        {   // Pick the first candidate that doesn't result in an edge that intersects any existing edges
            pCandidate = iter->second;
            auto start{pCandidate == pFirstHit ? std::next(concaveHull.begin()) : concaveHull.begin()};
            auto finish{std::prev(concaveHull.end())};

            intersects = false;
            for (auto hullIter = start; hullIter != finish && hullIter != concaveHull.end(); ++hullIter)
            {
                // Get the edge associated with the current hull vertex
                const CaloHit *h0{*hullIter};
                const CaloHit *h1{*std::next(hullIter)};
                if (h1 == pCurrentHit || h0 == pCurrentHit)
                    continue;
                intersects = LArClusterHelper::Intersects(h0->GetPositionVector(), h1->GetPositionVector(), pCurrentHit->GetPositionVector(),
                    pCandidate->GetPositionVector());

                if (intersects)
                    break;
            }
        }
        if (intersects)
        {   // No solution for current value of k, try again with a larger valuea
            concaveHull.clear();
            GetConcaveHull(pCluster, 2 * kk - 1, concaveHull);
            return;
        }
        // Update the baseline vector as the edge from the candidate hit to the previously added hit
        baseline = pCurrentHit->GetPositionVector() - pCandidate->GetPositionVector();
        pCurrentHit = pCandidate;
        concaveHull.push_back(pCurrentHit);
        const CartesianVector &hitPos{pCurrentHit->GetPositionVector()};

        if (xmin > hitPos.GetX())
            xmin = hitPos.GetX();
        if (zmin > hitPos.GetZ())
            zmin = hitPos.GetZ();
        if (xmax < hitPos.GetX())
            xmax = hitPos.GetX();
        if (zmax < hitPos.GetZ())
            zmax = hitPos.GetZ();

        // Remove the hit added to the hull from the list of hits to consider
        auto iter = std::find(inputCaloHitQueue.begin(), inputCaloHitQueue.end(), pCurrentHit);
        inputCaloHitQueue.erase(iter);
        ++step;
    }
    bool fullyContained{true};
    const CartesianVector externalPosMin{xmin - 1.f, 0.f, zmin - 1.f};
    const CartesianVector externalPosMax{xmax + 1.f, 0.f, zmax + 1.f};
    // Check that all hits reside within the hull - Jordan Curve Theorem: Define a ray starting outside the hull (simple bounding box is
    // sufficient to determine the start point) and terminating at a hit under consideration, an odd number of edge crossings implies the
    // hit is inside. Complication caused by frequent colinarity of hits, so need to check from different sides of the hull. If at least
    // one of these checks yields an odd number of crossings, the hit is inside.
    for (const CaloHit *pCaloHit : inputCaloHitQueue)
    {
        int crossingsMin{0};
        int crossingsMax{0};
        const CartesianVector &hitPos{pCaloHit->GetPositionVector()};
        auto start{concaveHull.begin()};
        auto finish{std::prev(concaveHull.end())};
        for (auto hullIter = start; hullIter != finish && hullIter != concaveHull.end(); ++hullIter)
        {
            // Get the edge associated with the current hull vertex
            const CaloHit *h0{*hullIter};
            const CaloHit *h1{*std::next(hullIter)};

            if (LArClusterHelper::Intersects(h0->GetPositionVector(), h1->GetPositionVector(), externalPosMin, hitPos, false))
                ++crossingsMin;
            if (LArClusterHelper::Intersects(h0->GetPositionVector(), h1->GetPositionVector(), externalPosMax, hitPos, false))
                ++crossingsMax;
        }
        // Get the edge connecting the first and last points in the hull to ensure the surface is closed
        if (!(crossingsMin % 2 || crossingsMax % 2))
        {
            fullyContained = false;
            break;
        }
    }
    if (!fullyContained && (kk + 1 < inputCaloHitList.size()))
    {   // Hull does not contain all hits, try again with a larger k
        concaveHull.clear();
        GetConcaveHull(pCluster, 2 * kk - 1, concaveHull);
    }
    if (concaveHull.front() == concaveHull.back())
        concaveHull.pop_back();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsCounterClockwise(const CartesianVector &a, const CartesianVector &b, const CartesianVector &c)
{
    return (c.GetZ()- a.GetZ())*(b.GetX() - a.GetX()) > (b.GetZ() - a.GetZ())*(c.GetX() - a.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::Intersects(const CartesianVector &a, const CartesianVector &b, const CartesianVector &c, const CartesianVector &d,
    const bool allowColinear)
{
    bool acdCCW{LArClusterHelper::IsCounterClockwise(a, c, d)};
    bool bcdCCW{LArClusterHelper::IsCounterClockwise(b, c, d)};
    bool abcCCW{LArClusterHelper::IsCounterClockwise(a, b, c)};
    bool abdCCW{LArClusterHelper::IsCounterClockwise(a, b, d)};
    bool diff1{acdCCW != bcdCCW};
    bool diff2{abcCCW != abdCCW};
    if (diff1 && diff2)
        return true;
    else if (diff1 || diff2)
        return allowColinear && LArClusterHelper::Colinear(a, b, d);
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::Colinear(const CartesianVector &a, const CartesianVector &b, const CartesianVector &c)
{
    const float dz1{b.GetZ() - c.GetZ()}, dz2{c.GetZ() - b.GetZ()}, dz3{a.GetZ() - b.GetZ()};
    const float area{a.GetX()*dz1 + b.GetX()*dz2 + c.GetX()*dz3};
    return std::fabs(area) < std::numeric_limits<float>::epsilon();
}

} // namespace lar_content
