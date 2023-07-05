/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/HitCorrelationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCorrelationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <chrono>
using namespace pandora;

namespace lar_content
{

HitCorrelationAlgorithm::HitCorrelationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HitCorrelationAlgorithm::~HitCorrelationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCorrelationAlgorithm::Run()
{
    m_volumeMap.clear();

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList, m_caloHitListName));

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
        if (!pLArCaloHit)
            continue;
        const unsigned int tpcVolume{pLArCaloHit->GetLArTPCVolumeId()};
        const unsigned int childVolume{pLArCaloHit->GetDaughterVolumeId()};
        const unsigned int key{1000 * tpcVolume + childVolume};
        m_volumeMap[key].AddCaloHit(pLArCaloHit);
        m_availabilityMap[pLArCaloHit] = true;
    }

    CaloHitList caloHitList3D;
    for (const auto & [ key, volume ] : m_volumeMap)
    {
        CaloHitList caloHitListU{volume.GetCaloHits(HitType::TPC_VIEW_U)};
        CaloHitList caloHitListV{volume.GetCaloHits(HitType::TPC_VIEW_V)};
        CaloHitList caloHitListW{volume.GetCaloHits(HitType::TPC_VIEW_W)};
        for (float scale : {1.f})//{0.1f, 0.25f, 0.5f, 1.0f})
        {
            std::cout << " Hits: " << caloHitListU.size() << " " << caloHitListV.size() << " " << caloHitListW.size() << std::endl;
            LArTripletVector triplets;
            auto start{std::chrono::high_resolution_clock::now()};
            this->Correlate(caloHitListU, caloHitListV, caloHitListW, scale, triplets);
            this->Correlate(caloHitListU, caloHitListW, caloHitListV, scale, triplets);
            this->Correlate(caloHitListV, caloHitListW, caloHitListU, scale, triplets);
            auto stop{std::chrono::high_resolution_clock::now()};
            std::cout << "Correlate " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;

            auto lambda{
                [](const LArTriplet &triplet1, const LArTriplet &triplet2) -> bool
                {
                    return std::get<3>(triplet1) < std::get<3>(triplet2);
                }
            };
            std::sort(triplets.begin(), triplets.end(), lambda);
            std::map<const CaloHit *, bool> usedHits;
            LArTripletVector filteredTriplets;
            for (const LArTriplet &triplet : triplets)
            {
                const CaloHit *pCaloHitU{std::get<0>(triplet)};
                const CaloHit *pCaloHitV{std::get<1>(triplet)};
                const CaloHit *pCaloHitW{std::get<2>(triplet)};
                if (usedHits.find(pCaloHitU) == usedHits.end() && usedHits.find(pCaloHitV) == usedHits.end() &&
                    usedHits.find(pCaloHitW) == usedHits.end())
                {
                    usedHits[pCaloHitU] = true;
                    usedHits[pCaloHitV] = true;
                    usedHits[pCaloHitW] = true;
                    filteredTriplets.emplace_back(triplet);
                }
            }
            std::cout << "Filtered triplet vector size: " << filteredTriplets.size() << std::endl;

            for (const LArTriplet &triplet : filteredTriplets)
            {
                const CaloHit *const pCaloHitU{std::get<0>(triplet)};
                const CaloHit *const pCaloHitV{std::get<1>(triplet)};
                const CaloHit *const pCaloHitW{std::get<2>(triplet)};
                const CartesianVector &posU{pCaloHitU->GetPositionVector()};
                const CartesianVector &posV{pCaloHitV->GetPositionVector()};
                const CartesianVector &posW{pCaloHitW->GetPositionVector()};
                CartesianVector pos3D(0, 0, 0);
                float chi2{0.f};
                LArGeometryHelper::MergeThreeWidePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, posU, posV, posW,
                    0.5f * pCaloHitU->GetCellSize1(), 0.5f * pCaloHitV->GetCellSize1(), 0.5f * pCaloHitW->GetCellSize1(), pos3D, chi2);
                const CaloHit *pCaloHit3D{nullptr};
                this->Create3DHit(triplet, pos3D, pCaloHit3D);

                if (!pCaloHit3D)
                    return STATUS_CODE_FAILURE;

                caloHitList3D.emplace_back(pCaloHit3D);
            }
        }
    }

    std::string caloHitListName{"CaloHitList3D"};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList3D, caloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCorrelationAlgorithm::Correlate(const CaloHitList &caloHitList1, const CaloHitList &caloHitList2, const CaloHitList &caloHitListTarget,
    const float scale, LArTripletVector &triplets) const
{
    const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
    if (caloHitList1.empty() || caloHitList2.empty() || caloHitListTarget.empty())
        return;

    const float minX{std::max(caloHitList1.front()->GetPositionVector().GetX(), caloHitList2.front()->GetPositionVector().GetX())};
    const float maxX{std::min(caloHitList1.back()->GetPositionVector().GetX(), caloHitList2.back()->GetPositionVector().GetX())};

    if (minX > maxX)
        return;

    const float chi2CutOff{1.f};
    for (auto iter = caloHitList1.begin(); iter != caloHitList1.end(); )
    {
        if ((*iter)->GetPositionVector().GetX() < minX)
        {
            ++iter;
            continue;
        }
        const CaloHit *const pCaloHit1{*iter};
        const float hitWidth1{0.5f * pCaloHit1->GetCellSize1()};
        const float x1{pCaloHit1->GetPositionVector().GetX()};
        if (x1 > maxX)
            break;
        const float lowX1{x1 - scale * hitWidth1};
        const float highX1{x1 + scale * hitWidth1};

        for (const CaloHit *pCaloHit2 : caloHitList2)
        {
            const float x2{pCaloHit2->GetPositionVector().GetX()};
            const float hitWidth2{0.5f * pCaloHit2->GetCellSize1()};
            const float lowX2{x2 - scale * hitWidth2};
            const float highX2{x2 + scale * hitWidth2};
            if (((x2 >= lowX1) && (x2 <= highX1)) || ((x1 >= lowX2) && (x1 <= highX2)))
            {
                const float lowXT{std::min(lowX1, lowX2)}, highXT{std::max(highX1, highX2)};
                if (pCaloHit1->GetHitType() == TPC_VIEW_U && pCaloHit2->GetHitType() == TPC_VIEW_V)
                {
                    const double u{pCaloHit1->GetPositionVector().GetZ()};
                    const double v{pCaloHit2->GetPositionVector().GetZ()};
                    const double w{transform->UVtoW(u, v)};
                    // Find the closest hit in W
                    double minDistance{std::numeric_limits<float>::max()};
                    const CaloHit *pClosestHit{nullptr};
                    for (const CaloHit *pCaloHitT : caloHitListTarget)
                    {
                        const float xT{pCaloHitT->GetPositionVector().GetX()};
                        if ((xT >= lowXT) && (xT <= highXT))
                        {
                            const float wT{pCaloHitT->GetPositionVector().GetZ()};
                            const double distance2{(xT - 0.5 * (x1 + x2)) * (xT - 0.5 * (x1 + x2)) + (wT - w) * (wT - w)};
                            if (distance2 < minDistance)
                            {
                                minDistance = distance2;
                                pClosestHit = pCaloHitT;
                            }
                        }
                    }
                    if (pClosestHit)
                    {
                        const double sigmaU{0.5 * pCaloHit1->GetCellSize0()};
                        const double sigmaV{0.5 * pCaloHit2->GetCellSize0()};
                        const double sigmaW{0.5 * pClosestHit->GetCellSize0()};
                        const float wT{pClosestHit->GetPositionVector().GetZ()};
                        double y{0.}, z{0.}, chi2{0.};
                        transform->GetMinChiSquaredYZ(u, v, wT, sigmaU, sigmaV, sigmaW, y, z, chi2);
                        if (chi2 < chi2CutOff)
                        {
                            const LArTriplet hitTriplet{std::make_tuple(pCaloHit1, pCaloHit2, pClosestHit, chi2)};
                            triplets.emplace_back(hitTriplet);
                        }
                    }
                }
                else if (pCaloHit1->GetHitType() == TPC_VIEW_U && pCaloHit2->GetHitType() == TPC_VIEW_W)
                {
                    const double u{pCaloHit1->GetPositionVector().GetZ()};
                    const double w{pCaloHit2->GetPositionVector().GetZ()};
                    const double v{transform->WUtoV(w, u)};
                    // Find the closest hit in W
                    float minDistance{std::numeric_limits<float>::max()};
                    const CaloHit *pClosestHit{nullptr};
                    for (const CaloHit *pCaloHitT : caloHitListTarget)
                    {
                        const float xT{pCaloHitT->GetPositionVector().GetX()};
                        if ((xT >= lowXT) && (xT <= highXT))
                        {
                            const float vT{pCaloHitT->GetPositionVector().GetZ()};
                            const double distance2{(xT - 0.5 * (x1 + x2)) * (xT - 0.5 * (x1 + x2)) + (vT - v) * (vT - v)};
                            if (distance2 < minDistance)
                            {
                                minDistance = distance2;
                                pClosestHit = pCaloHitT;
                            }
                        }
                    }
                    if (pClosestHit)
                    {
                        const double sigmaU{0.5 * pCaloHit1->GetCellSize0()};
                        const double sigmaW{0.5 * pCaloHit2->GetCellSize0()};
                        const double sigmaV{0.5 * pClosestHit->GetCellSize0()};
                        const float vT{pClosestHit->GetPositionVector().GetZ()};
                        double y{0.}, z{0.}, chi2{0.};
                        transform->GetMinChiSquaredYZ(u, vT, w, sigmaU, sigmaV, sigmaW, y, z, chi2);
                        if (chi2 < chi2CutOff)
                        {
                            const LArTriplet hitTriplet{std::make_tuple(pCaloHit1, pClosestHit, pCaloHit2, chi2)};
                            triplets.emplace_back(hitTriplet);
                        }
                    }
                }
                if (pCaloHit1->GetHitType() == TPC_VIEW_V && pCaloHit2->GetHitType() == TPC_VIEW_W)
                {
                    const double v{pCaloHit1->GetPositionVector().GetZ()};
                    const double w{pCaloHit2->GetPositionVector().GetZ()};
                    const double u{transform->VWtoU(v, w)};
                    // Find the closest hit in W
                    float minDistance{std::numeric_limits<float>::max()};
                    const CaloHit *pClosestHit{nullptr};
                    for (const CaloHit *pCaloHitT : caloHitListTarget)
                    {
                        const float xT{pCaloHitT->GetPositionVector().GetX()};
                        if ((xT >= lowXT) && (xT <= highXT))
                        {
                            const float uT{pCaloHitT->GetPositionVector().GetZ()};
                            const double distance2{(xT - 0.5 * (x1 + x2)) * (xT - 0.5 * (x1 + x2)) + (uT - u) * (uT - u)};
                            if (distance2 < minDistance)
                            {
                                minDistance = distance2;
                                pClosestHit = pCaloHitT;
                            }
                        }
                    }
                    if (pClosestHit)
                    {
                        const double sigmaV{0.5 * pCaloHit1->GetCellSize0()};
                        const double sigmaW{0.5 * pCaloHit2->GetCellSize0()};
                        const double sigmaU{0.5 * pClosestHit->GetCellSize0()};
                        const float uT{pClosestHit->GetPositionVector().GetZ()};
                        double y{0.}, z{0.}, chi2{0.};
                        transform->GetMinChiSquaredYZ(uT, v, w, sigmaU, sigmaV, sigmaW, y, z, chi2);
                        if (chi2 < chi2CutOff)
                        {
                            const LArTriplet hitTriplet{std::make_tuple(pClosestHit, pCaloHit1, pCaloHit2, chi2)};
                            triplets.emplace_back(hitTriplet);
                        }
                    }
                }
            }
        }
        ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCorrelationAlgorithm::Create3DHit(const LArTriplet &hitTriplet, const CartesianVector &pos, const CaloHit *&pCaloHit3D) const
{
    LArCaloHit *pCaloHitU{dynamic_cast<LArCaloHit *>(const_cast<CaloHit *>(std::get<0>(hitTriplet)))};
    LArCaloHit *pCaloHitV{dynamic_cast<LArCaloHit *>(const_cast<CaloHit *>(std::get<1>(hitTriplet)))};
    LArCaloHit *pCaloHitW{dynamic_cast<LArCaloHit *>(const_cast<CaloHit *>(std::get<2>(hitTriplet)))};
    if (pCaloHitU && pCaloHitV && pCaloHitW)
    {
        pCaloHitU->AddSiblingHit(pCaloHitV);
        pCaloHitU->AddSiblingHit(pCaloHitW);
        pCaloHitV->AddSiblingHit(pCaloHitU);
        pCaloHitV->AddSiblingHit(pCaloHitW);
        pCaloHitW->AddSiblingHit(pCaloHitU);
        pCaloHitW->AddSiblingHit(pCaloHitV);
    }

    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pos;
    parameters.m_hitType = TPC_3D;

    const CaloHit *const pCaloHit2D(std::get<2>(hitTriplet));
    parameters.m_pParentAddress = static_cast<const void *>(pCaloHit2D);

    parameters.m_cellThickness = pCaloHit2D->GetCellThickness();
    parameters.m_cellGeometry = RECTANGULAR;
    parameters.m_cellSize0 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellSize1 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellNormalVector = pCaloHit2D->GetCellNormalVector();
    parameters.m_expectedDirection = pCaloHit2D->GetExpectedDirection();
    parameters.m_nCellRadiationLengths = pCaloHit2D->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit2D->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit2D->GetTime();
    parameters.m_inputEnergy = pCaloHit2D->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit2D->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit2D->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit2D->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit2D->IsDigital();
    parameters.m_hitRegion = pCaloHit2D->GetHitRegion();
    parameters.m_layer = pCaloHit2D->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit2D->IsInOuterSamplingLayer();
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCorrelationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void HitCorrelationAlgorithm::TpcChildVolume::AddCaloHit(const CaloHit *pCaloHit)
{
    switch (pCaloHit->GetHitType())
    {
        case HitType::TPC_VIEW_U:
            m_hitsU.emplace_back(pCaloHit);
            break;
        case HitType::TPC_VIEW_V:
            m_hitsV.emplace_back(pCaloHit);
            break;
        default:
            m_hitsW.emplace_back(pCaloHit);
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList HitCorrelationAlgorithm::TpcChildVolume::GetCaloHits(const HitType type) const
{
    switch (type)
    {
        case HitType::TPC_VIEW_U:
        {
            std::vector<const CaloHit *> sortedHitList;
            for (const CaloHit *pCaloHit : m_hitsU)
                sortedHitList.emplace_back(pCaloHit);
            std::sort(sortedHitList.begin(), sortedHitList.end(), LArClusterHelper::SortHitsByPositionInX);
            CaloHitList caloHitList;
            for (const CaloHit *pCaloHit : sortedHitList)
                caloHitList.emplace_back(pCaloHit);
            return caloHitList;
        }
        case HitType::TPC_VIEW_V:
        {
            std::vector<const CaloHit *> sortedHitList;
            for (const CaloHit *pCaloHit : m_hitsV)
                sortedHitList.emplace_back(pCaloHit);
            std::sort(sortedHitList.begin(), sortedHitList.end(), LArClusterHelper::SortHitsByPositionInX);
            CaloHitList caloHitList;
            for (const CaloHit *pCaloHit : sortedHitList)
                caloHitList.emplace_back(pCaloHit);
            return caloHitList;
        }
        default:
        {
            std::vector<const CaloHit *> sortedHitList;
            for (const CaloHit *pCaloHit : m_hitsW)
                sortedHitList.emplace_back(pCaloHit);
            std::sort(sortedHitList.begin(), sortedHitList.end(), LArClusterHelper::SortHitsByPositionInX);
            CaloHitList caloHitList;
            for (const CaloHit *pCaloHit : sortedHitList)
                caloHitList.emplace_back(pCaloHit);
            return caloHitList;
        }
    }
}

} // namespace lar_content
