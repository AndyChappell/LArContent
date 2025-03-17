/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include <Eigen/Dense>

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/KalmanClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

KalmanClusterCreationAlgorithm::KalmanClusterCreationAlgorithm() :
    m_minMipFraction{0.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanClusterCreationAlgorithm::~KalmanClusterCreationAlgorithm()
{
    this->Reset();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::Run()
{
    this->Reset();
    std::map<HitType, int> viewHitCountMap;
    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        m_viewHitsMap[view] = OrderedCaloHitList();
        viewHitCountMap[view] = 0;
    }
    // Collect hits from all views into layer-ordered hit lists
    for (std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (!pCaloHitList->empty())
        {
            HitType view{pCaloHitList->front()->GetHitType()};
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FilterCaloHits(pCaloHitList, m_viewHitsMap[view]));
            for (OrderedCaloHitList::const_iterator iter = m_viewHitsMap[view].begin(), iterEnd = m_viewHitsMap[view].end(); iter != iterEnd; ++iter)
                viewHitCountMap[view] += iter->second->size();
        }
    }

    float min{std::numeric_limits<float>::max()}, max{std::numeric_limits<float>::lowest()};
    this->GetSpanX(min, max);
    for (std::string listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (!pCaloHitList->empty())
        {
            HitType view{pCaloHitList->front()->GetHitType()};
            if (m_slicedCaloHits.find(view) == m_slicedCaloHits.end())
            {
                m_slicedCaloHits[view] = new LArSlicedCaloHitList(*pCaloHitList, min, max);
            }
            else
            {
                std::cout << "Error: KalmanClusterCreationAlgorithm - Duplciate view specified in calo hit lists" << std::endl;
                throw StatusCodeException(STATUS_CODE_ALREADY_INITIALIZED);
            }
        }
    }
    for (HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        if (m_slicedCaloHits.find(view) == m_slicedCaloHits.end())
            m_slicedCaloHits[view] = new LArSlicedCaloHitList(CaloHitList(), min, max);
    }

    // Order the views by hit count, because we want to use the view with the most hits as our "primary" view
    std::vector<HitType> order;
    if (viewHitCountMap[TPC_VIEW_W] >= viewHitCountMap[TPC_VIEW_V])
    {
        order.emplace_back(TPC_VIEW_W);
        order.emplace_back(TPC_VIEW_V);
        if (viewHitCountMap[TPC_VIEW_U] > viewHitCountMap[TPC_VIEW_W])
            order.insert(order.begin(), TPC_VIEW_U);
        else if (viewHitCountMap[TPC_VIEW_U] <= viewHitCountMap[TPC_VIEW_V])
            order.insert(order.end(), TPC_VIEW_U);
        else
            order.insert(order.begin() + 1, TPC_VIEW_U);
    }
    else
    {
        order.emplace_back(TPC_VIEW_V);
        order.emplace_back(TPC_VIEW_W);
        if (viewHitCountMap[TPC_VIEW_U] > viewHitCountMap[TPC_VIEW_V])
            order.insert(order.begin(), TPC_VIEW_U);
        else if (viewHitCountMap[TPC_VIEW_U] <= viewHitCountMap[TPC_VIEW_W])
            order.insert(order.end(), TPC_VIEW_U);
        else
            order.insert(order.begin() + 1, TPC_VIEW_U);
    }

    for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        for (OrderedCaloHitList::const_iterator iter = m_viewHitsMap[view].begin(), iterEnd = m_viewHitsMap[view].end(); iter != iterEnd; ++iter)
        {
            CaloHitVector caloHits(iter->second->begin(), iter->second->end());
            std::sort(caloHits.begin(), caloHits.end(), LArClusterHelper::SortHitsByPosition);
        }
    }
    if (viewHitCountMap[order.back()] > 0)
    {
        // All three views have hits
        this->IdentifyCandidateClusters(order);
    }
    else if (viewHitCountMap[order[1]] > 0)
    {
        // Only two views have hits - suspect we won't need this special case
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::IdentifyCandidateClusters(const ViewVector &order)
{
    // Begin processing the "primary" view with reference to the secondary and tertiary views
    const LArSlicedCaloHitList::SliceHitMap &sliceHitMap0{m_slicedCaloHits[order.at(0)]->GetSliceHitMap()};
    const LArSlicedCaloHitList::SliceHitMap &sliceHitMap1{m_slicedCaloHits[order.at(1)]->GetSliceHitMap()};
    const LArSlicedCaloHitList::SliceHitMap &sliceHitMap2{m_slicedCaloHits[order.at(2)]->GetSliceHitMap()};
    for (const auto &[bin, caloHits0] : sliceHitMap0)
    {
        CaloHitVector caloHits1, caloHits2;
        this->GetSlices(sliceHitMap1, bin, caloHits1);
        this->GetSlices(sliceHitMap2, bin, caloHits2);

        // Create a vector of all hit permutations involing a hit from caloHits0
        const size_t n3HitTuples{caloHits0.size() * caloHits1.size() * caloHits2.size()};
        const size_t n2HitTuples{caloHits0.size() * caloHits1.size() + caloHits0.size() * caloHits2.size()};
        CandidateCluster::HitTripletVector triplets(n3HitTuples + n2HitTuples, std::make_tuple(nullptr, nullptr, nullptr));
        CartesianPointVector hits3D;
        FloatVector chi2s;
        this->Make3DHitPermutations(caloHits0, caloHits1, caloHits2, triplets);
        this->Filter3DHitPermutations(triplets, hits3D, chi2s);
        std::cout << "Number of 3D hit permutations: " << hits3D.size() << std::endl;

        // Now define the full set of possible two hit starting clusters and try to build on those and see which is best
        typedef std::tuple<size_t, size_t, CartesianVector, CartesianVector> SeedPair;
        std::vector<SeedPair> seedPairs;
        for (size_t i = 0; i < hits3D.size(); ++i)
        {
            // Only make seeds from triplets that are fully populated
            const CaloHit *const pCaloHit0{std::get<2>(triplets.at(i))};
            if (!pCaloHit0)
                continue;
            for (size_t j = i + 1; j < hits3D.size(); ++j)
            {
                const CaloHit *const pCaloHit1{std::get<2>(triplets.at(j))};
                if (!pCaloHit1)
                    continue;
                seedPairs.emplace_back(std::make_tuple(i, j, hits3D.at(i), hits3D.at(j)));
            }
        }

        std::cout << "Number of seed pairs: " << seedPairs.size() << std::endl;

        // Now we have a list of 3D hit positions and chi2 values, so we can start clustering
        for (auto seed : seedPairs)
        {
            const size_t s1{std::get<0>(seed)};
            const size_t s2{std::get<1>(seed)};
            const CartesianVector &seed1{std::get<2>(seed)};
            Eigen::VectorXd init(3);
            init << seed1.GetX(), seed1.GetY(), seed1.GetZ();
            KalmanFilter3D kalmanFilter(0.5, 0.1, 1.0, init);
            kalmanFilter.Predict();
            const CartesianVector &seed2{std::get<3>(seed)};
            Eigen::VectorXd init2(3);
            init2 << seed2.GetX(), seed2.GetY(), seed2.GetZ();
            kalmanFilter.Update(init2);
            std::cout << "Kalman filter initialized with positions: (" << init.coeff(0) << " " << init.coeff(1) << " " << init.coeff(2) << ") (" <<
                init2.coeff(0) << " " << init2.coeff(1) << " " << init2.coeff(2) << ")" << std::endl;

            for (size_t i = 0; i < hits3D.size(); ++i)
            {
                if (i == s1 || i == s2)
                    continue;
                const CartesianVector &other{hits3D.at(i)};
                Eigen::VectorXd measurement(3);
                measurement << other.GetX(), other.GetY(), other.GetZ();
                std::cout << "   Kalman filter measurement: (" << measurement.coeff(0) << " " << measurement.coeff(1) << " " << measurement.coeff(2) << ")" << std::endl;
                kalmanFilter.Predict();
                const Eigen::VectorXd &state{kalmanFilter.GetTemporaryState()};
                std::cout << "   Kalman filter prediction: (" << state.coeff(0) << " " << state.coeff(1) << " " << state.coeff(2) << ")" << std::endl;
                //kalmanFilter.Update(measurement);
                /*
                const float chi2{chi2s.at(i)};
                const CandidateCluster::HitTriplet &triplet{triplets.at(i)};
                const CaloHit *const pHit0{std::get<0>(triplet)};
                const CaloHit *const pHit1{std::get<1>(triplet)};
                const CaloHit *const pHit2{std::get<2>(triplet)};
                */
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::Make3DHitPermutations(const CaloHitVector &caloHits0, const CaloHitVector &caloHits1, const CaloHitVector &caloHits2, CandidateCluster::HitTripletVector &triplets) const
{
    int i{0};
    // Construct 3 hit tuples
    for (const CaloHit *pCaloHit0 : caloHits0)
    {
        for (const CaloHit *pCaloHit1 : caloHits1)
        {
            for (const CaloHit *pCaloHit2 : caloHits2)
            {
                std::get<0>(triplets[i]) = pCaloHit0;
                std::get<1>(triplets[i]) = pCaloHit1;
                std::get<2>(triplets[i]) = pCaloHit2;
                ++i;
            }
        }
    }
    // Construct 2 hit tuples
    for (const CaloHit *pCaloHit0 : caloHits0)
    {
        for (const CaloHit *pCaloHit1 : caloHits1)
        {
            std::get<0>(triplets[i]) = pCaloHit0;
            std::get<1>(triplets[i]) = pCaloHit1;
            ++i;
        }
    }
    for (const CaloHit *pCaloHit0 : caloHits0)
    {
        for (const CaloHit *pCaloHit2 : caloHits2)
        {
            std::get<0>(triplets[i]) = pCaloHit0;
            std::get<1>(triplets[i]) = pCaloHit2;
            ++i;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::Filter3DHitPermutations(CandidateCluster::HitTripletVector &triplets, CartesianPointVector &hits3D, FloatVector &chi2s) const
{
    for (auto iter = triplets.begin(); iter != triplets.end();)
    {
        const CandidateCluster::HitTriplet &triplet{*iter};
        const CaloHit *const pHit0{std::get<0>(triplet)};
        const CaloHit *const pHit1{std::get<1>(triplet)};
        const CaloHit *const pHit2{std::get<2>(triplet)};
        CartesianVector pos3D(0, 0, 0);
        float chi2{0.f};
        if (pHit2)
        {
            LArGeometryHelper::MergeThreeWideHits3D(this->GetPandora(), *pHit0, *pHit1, *pHit2, pos3D, chi2);
            if (chi2 <= 6.f)
            {
                ++iter;
                hits3D.emplace_back(pos3D);
                chi2s.emplace_back(chi2);
            }
            else
            {
                iter = triplets.erase(iter);
            }
        }
        else
        {
            LArGeometryHelper::MergeTwoWideHits3D(this->GetPandora(), *pHit0, *pHit1, pos3D, chi2);
            if (chi2 <= 2.f)
            {
                ++iter;
                hits3D.emplace_back(pos3D);
                chi2s.emplace_back(chi2);
            }
            else
            {
                iter = triplets.erase(iter);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::FilterCaloHits(const CaloHitList *const pCaloHitList, OrderedCaloHitList &selectedCaloHitList) const
{
    CaloHitList availableHitList;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (PandoraContentApi::IsAvailable(*this, pCaloHit) && pCaloHit->GetMipEquivalentEnergy() >= m_minMipFraction)
            availableHitList.push_back(pCaloHit);
    }

    if (availableHitList.empty())
        return STATUS_CODE_SUCCESS;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, selectedCaloHitList.Add(availableHitList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::GetSpanX(float &min, float &max) const
{
    min = std::numeric_limits<float>::max();
    max = std::numeric_limits<float>::lowest();

    for (const auto &[view, orderedHitList] : m_viewHitsMap)
    {
        for (OrderedCaloHitList::const_iterator iter = orderedHitList.begin(), iterEnd = orderedHitList.end(); iter != iterEnd; ++iter)
        {
            const CaloHitList &caloHitList{*iter->second};
            for (const CaloHit *const pCaloHit : caloHitList)
            {
                const float xLo{pCaloHit->GetPositionVector().GetX() - 0.5f * pCaloHit->GetCellSize1()};
                const float xHi{pCaloHit->GetPositionVector().GetX() + 0.5f * pCaloHit->GetCellSize1()};
                if (xLo < min)
                    min = xLo;
                if (xHi > max)
                    max = xHi;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::Reset()
{
    StatusCode status{Algorithm::Reset()};
    for (auto &[bin, slicedCaloHitList] : m_slicedCaloHits)
    {
        delete slicedCaloHitList;
    }
    m_slicedCaloHits.clear();
    m_viewHitsMap.clear();

    return status;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::GetSlices(const LArSlicedCaloHitList::SliceHitMap &sliceHitMap, const size_t bin, CaloHitVector &caloHits_m1,
    CaloHitVector &caloHits_p1) const
{
    if (sliceHitMap.find(bin - 1) != sliceHitMap.end())
    {
        const CaloHitVector &caloHits{sliceHitMap.at(bin - 1)};
        caloHits_m1.insert(caloHits_m1.end(), caloHits.begin(), caloHits.end());
    }
    if (sliceHitMap.find(bin + 1) != sliceHitMap.end())
    {
        const CaloHitVector &caloHits{sliceHitMap.at(bin + 1)};
        caloHits_p1.insert(caloHits_p1.end(), caloHits.begin(), caloHits.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::GetSlices(const LArSlicedCaloHitList::SliceHitMap &sliceHitMap, const size_t bin, CaloHitVector &caloHits) const
{
    for (int i = -1; i <= 1; ++i)
    {
        if (sliceHitMap.find(bin + i) != sliceHitMap.end())
        {
            const CaloHitVector &sliceCaloHits{sliceHitMap.at(bin + i)};
            caloHits.insert(caloHits.end(), sliceCaloHits.begin(), sliceCaloHits.end());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMipFraction", m_minMipFraction));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet(const CaloHit *const pHit1, const CaloHit *const pHit2, const CaloHit *const pHit3)
{
    const CaloHit *pHitU{nullptr}, *pHitV{nullptr}, *pHitW{nullptr};
    if (pHit1)
    {
        switch (pHit1->GetHitType())
        {
            case TPC_VIEW_U:
                pHitU = pHit1;
                break;
            case TPC_VIEW_V:
                pHitU = pHit1;
                break;
            case TPC_VIEW_W:
                pHitU = pHit1;
                break;
            default:
                std::cout << "Error: KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet - Unexpetced hit type: " << pHit1->GetHitType() << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }
    if (pHit2)
    {
        switch (pHit2->GetHitType())
        {
            case TPC_VIEW_U:
                pHitU = pHit2;
                break;
            case TPC_VIEW_V:
                pHitU = pHit2;
                break;
            case TPC_VIEW_W:
                pHitU = pHit2;
                break;
            default:
                std::cout << "Error: KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet - Unexpetced hit type: " << pHit2->GetHitType() << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }
    if (pHit3)
    {
        switch (pHit3->GetHitType())
        {
            case TPC_VIEW_U:
                pHitU = pHit3;
                break;
            case TPC_VIEW_V:
                pHitU = pHit3;
                break;
            case TPC_VIEW_W:
                pHitU = pHit3;
                break;
            default:
                std::cout << "Error: KalmanClusterCreationAlgorithm::CandidateCluster::AddTriplet - Unexpetced hit type: " << pHit3->GetHitType() << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }
    m_hitTriplets.emplace_back(std::make_tuple(pHitU, pHitV, pHitW));
}

} // namespace lar_content
