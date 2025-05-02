/**
 *  @file   larpandoracontent/LArHelpers/LArSlicedCaloHitList.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArObjects/LArSlicedCaloHitList.h"

using namespace pandora;

namespace lar_content
{

LArSlicedCaloHitList::LArSlicedCaloHitList(const CaloHitList &caloHitList, const float xMin, const float xMax, const float sliceSize) :
    m_xMin{xMin},
    m_xMax{xMax + std::numeric_limits<float>::epsilon()},
    m_sliceSize{sliceSize}
{
    if (m_sliceSize < std::numeric_limits<float>::epsilon())
    {
        m_sliceSize = 1.f;
        std::cerr << "WARNING: LArSlicedCaloHitList slice size less than " << std::numeric_limits<float>::epsilon() << ". Defaulting to 1.0f" << std::endl;
    }
    m_nBins = static_cast<size_t>(std::ceil((m_xMax - m_xMin) / m_sliceSize));

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float halfWidth{0.5f * pCaloHit->GetCellSize1()};
        const float lo{x - halfWidth};
        const float hi{x + halfWidth};
        for (float pos = lo; pos <= hi; pos += m_sliceSize)
        {
            const int bin{static_cast<int>((pos - m_xMin) / m_sliceSize)};
            m_sliceHitMap[bin].emplace_back(pCaloHit);
        }
    }
    for (auto &[bin, caloHits] : m_sliceHitMap)
        std::sort(caloHits.begin(), caloHits.end(), LArClusterHelper::SortHitsByPositionInX);

    // Ensure unpopulated bins within the [min, max] region are allocated an empty hit vector
    for (float x = m_xMin; x < m_xMax; x += m_sliceSize)
    {
        const int bin{static_cast<int>((x - m_xMin) / m_sliceSize)};
        if (m_sliceHitMap.find(bin) == m_sliceHitMap.end())
            m_sliceHitMap[bin] = CaloHitVector();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArSlicedCaloHitList::~LArSlicedCaloHitList()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitVector &LArSlicedCaloHitList::GetCaloHitsInSlice(const size_t slice) const
{
    if (this->find(slice) != m_sliceHitMap.end())
        return m_sliceHitMap.at(slice);
    else
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

size_t LArSlicedCaloHitList::GetNCaloHitsInSlice(const size_t slice) const
{
    if (this->find(slice) != m_sliceHitMap.end())
        return m_sliceHitMap.at(slice).size();
    else
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArSlicedCaloHitList::SliceHitMap &LArSlicedCaloHitList::GetSliceHitMap() const
{
    return m_sliceHitMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

size_t LArSlicedCaloHitList::GetSliceFromCoordinate(const float x) const
{
    if (x < m_xMin || x > m_xMax)
    {
        std::cout << "Error: LArSlicedCaloHitList - Requested slice for x = " << x << " Range: [" << m_xMin << ","  << m_xMax << "]" << std::endl;
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
    }

    return static_cast<size_t>((x - m_xMin) / m_sliceSize);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArSlicedCaloHitList::const_iterator LArSlicedCaloHitList::begin() const
{
    return m_sliceHitMap.begin();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArSlicedCaloHitList::const_iterator LArSlicedCaloHitList::end() const
{
    return m_sliceHitMap.end();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArSlicedCaloHitList::const_reverse_iterator LArSlicedCaloHitList::rbegin() const
{
    return m_sliceHitMap.rbegin();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArSlicedCaloHitList::const_reverse_iterator LArSlicedCaloHitList::rend() const
{
    return m_sliceHitMap.rend();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArSlicedCaloHitList::const_iterator LArSlicedCaloHitList::find(const size_t slice) const
{
    return m_sliceHitMap.find(slice);
}

//------------------------------------------------------------------------------------------------------------------------------------------

size_t LArSlicedCaloHitList::size() const
{
    return m_sliceHitMap.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArSlicedCaloHitList::empty() const
{
    return m_sliceHitMap.empty();
}

} // namespace lar_content
