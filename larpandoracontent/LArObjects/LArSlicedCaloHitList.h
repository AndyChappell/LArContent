/**
 *  @file   larpandoracontent/LArHelpers/LArSlicedCaloHitList.h
 *
 *  @brief  Header file for the LArSlicedCaloHitList class.
 *
 *  $Log: $
 */
#ifndef LAR_SLICED_CALO_HIT_LIST_H
#define LAR_SLICED_CALO_HIT_LIST_H 1

#include "Objects/CaloHit.h"

#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include <map>

namespace lar_content
{

/**
 *  @brief  LArSlicedCaloHitList class
 */
class LArSlicedCaloHitList
{
public:
    typedef std::map<size_t, pandora::CaloHitVector> SliceHitMap;
    typedef SliceHitMap::value_type value_type;
    typedef SliceHitMap::const_iterator const_iterator;
    typedef SliceHitMap::const_reverse_iterator const_reverse_iterator;

    /**
     *  @brief  Constructs a representation of a calo hit list where hits are partitioned into slices in x.
     *          ATTN: This is intended as a convenient lookup structure, and so the pointer to a hit with a width that spans multiple slices
     *          will be stored in each of those slice lists.
     *
     *          The slices are constructed based on the minimum and maximum x coordinates and a (default) 1 cm slice width. To maximise the
     *          usefulness/efficiency of this class, these minima and maxima should be defined consistently across views, rather than defined
     *          within a single view.
     *
     *  @param  caloHitList fullyConnect Whether or not to connect disconnected regions
     *  @param  xMin the minimum x coordinate to consider
     *  @param  xMax the maximum x coordinate to consider
     *  @param  sliceSize the size of a slice
     */
    LArSlicedCaloHitList(const pandora::CaloHitList &caloHitList, const float xMin, const float xMax, const float sliceSize = 1.f);

    /**
     *  @brief  Destructor
     */
    ~LArSlicedCaloHitList();

    /**
     *  @brief  Get the calo hits in a specified slice
     * 
     *  @param  slice the slice
     * 
     *  @return The calo hit vector in the slice
     */
    const pandora::CaloHitVector &GetCaloHitsInSlice(const size_t slice) const;

    /**
     *  @brief  Get the number of calo hits in a specified slice
     * 
     *  @param  slice the slice
     * 
     *  @return The number of calo hits in the specified slice
     */
    size_t GetNCaloHitsInSlice(const size_t slice) const;

    /**
     *  @brief  Get the map between slices and hits
     * 
     *  @return A const reference to the map between slices and hits
     */
    const SliceHitMap &GetSliceHitMap() const;

    /**
     *  @brief  Get the slice corresponding to a drift coordinate. Throws an exception if the coordinate is out of bounds.
     * 
     *  @param  x the coordinate whose slice is to be determined
     *
     *  @return The slice corresponding to the given drift coordinate
     */
    size_t GetSliceFromCoordinate(const float x) const;

    /**
     *  @brief  Returns a const iterator referring to the first element in the sliced calo hit list
     */
    const_iterator begin() const;

    /**
     *  @brief  Returns a const iterator referring to the past-the-end element in the sliced calo hit list
     */
    const_iterator end() const;

    /**
     *  @brief  Returns a const reverse iterator referring to the first element in the sliced calo hit list
     */
    const_reverse_iterator rbegin() const;

    /**
     *  @brief  Returns a const reverse iterator referring to the past-the-end element in the sliced calo hit list
     */
    const_reverse_iterator rend() const;

    /**
     *  @brief  Searches the container for an element with specified slice and returns an iterator to it if found,
     *          otherwise it returns an iterator to the past-the-end element.
     *
     *  @param  slice the slice to find
     */
    const_iterator find(const size_t slice) const;

    /**
     *  @brief  Returns the number of elements in the container.
     */
    size_t size() const;

    /**
     *  @brief  Returns whether the map container is empty
     */
    bool empty() const;

private:
    SliceHitMap m_sliceHitMap;  ///< The edges defining the graph
    float m_xMin;   ///< The minimum x-coordinate (including hit widths) for the slice
    float m_xMax;   ///< The maximum x-coordinate (including hit widths) for the slice
    float m_sliceSize;  ///< The size of each slice
    size_t m_nBins; ///< The number of bins in the map
};

} // namespace lar_content

#endif // #ifndef LAR_SLICED_CALO_HIT_LIST_H
