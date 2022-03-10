/**
 *  @file   larpandoracontent/LArUtility/TpcHitVolume.cc
 *
 *  @brief  Implementation of the list merging algorithm class.
 *
 *  $Log: $
 */

#include "Objects/Cluster.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArUtility/TpcHitVolume.h"

using namespace pandora;

namespace lar_content
{

TpcHitVolume::TpcHitVolume(const unsigned int cryostat, const unsigned int tpc, const unsigned int child, const pandora::CartesianVector &center,
    const pandora::CartesianVector &length, const LArTransformationPlugin *const pTransform) : m_cryostat{cryostat}, m_tpc{tpc}, m_child{child},
    m_min{center - length * 0.5}, m_max{center + length * 0.5}, m_uMin{0., 0., 0.}, m_uMax{0., 0., 0.}, m_vMin{0., 0., 0.}, m_vMax{0., 0., 0.}
{
    // ATTN: Wire planes with non-zero wire angles can map hits to virtual wire centres that extend beyond the physical bounds of the
    // daughter volume in the Z direction. As a result, minima and maximum for these planes must be adjusted to reflect this.
    const double uMinZ{pTransform->YZtoU(m_max.GetY(), m_min.GetZ())};
    const double uMaxZ{pTransform->YZtoU(m_min.GetY(), m_max.GetZ())};

    m_uMin = CartesianVector(m_min.GetX(), m_min.GetY(), uMinZ);
    m_uMax = CartesianVector(m_max.GetX(), m_max.GetY(), uMaxZ);

    const double vMinZ{pTransform->YZtoV(m_min.GetY(), m_min.GetZ())};
    const double vMaxZ{pTransform->YZtoV(m_max.GetY(), m_max.GetZ())};
    m_vMin = CartesianVector(m_min.GetX(), m_min.GetY(), vMinZ);
    m_vMax = CartesianVector(m_max.GetX(), m_max.GetY(), vMaxZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TpcHitVolume::~TpcHitVolume()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::Add(const CaloHitList &caloHitList)
{
    for (const CaloHit *const pCaloHit : caloHitList)
        this->Add(pCaloHit);
/*
    CartesianPointVector localCoords;
    if (!caloHitList.empty())
        this->GetLocalCoordinates(TPC_VIEW_V, localCoords);

    const CaloHitList caloHitListV{m_viewToCaloHitMap[TPC_VIEW_V]};
    auto hitIter{caloHitListV.begin()};
    auto posIter{localCoords.begin()};
    for (size_t i = 0; i < caloHitListV.size(); ++i)
    {
        const CartesianVector &local{*posIter};
        const CartesianVector &global{(*hitIter)->GetPositionVector()};
        std::cout << "TPC " << m_tpc << ":" << m_child << " " << (*hitIter)->GetHitType() << " Global: (" <<
            global.GetX() << "," << global.GetZ() << ") Local: (" << local.GetX() << "," << local.GetZ() << ")" << std::endl;
        ++posIter;
        ++hitIter;
    }*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::Add(const CaloHit *const pCaloHit)
{
    if (this->Contains(pCaloHit))
    {
        const HitType view{pCaloHit->GetHitType()};
        this->InitialiseView(view);
        m_viewToCaloHitMap[view].emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::GetLocalCoordinates(const HitType view, CartesianPointVector &localCoords) const
{
    // Probably want to check initialisation status here
    if (m_viewToCaloHitMap.find(view) != m_viewToCaloHitMap.end())
    {
        const CaloHitList &caloHitList{m_viewToCaloHitMap.at(view)};
        for (const CaloHit *pCaloHit : caloHitList)
            this->GetLocalCoordinate(pCaloHit, localCoords);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::GetLocalCoordinates(const Cluster *const pCluster, CartesianPointVector &localCoords) const
{
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    const CaloHitList &isolatedHitList{pCluster->GetIsolatedCaloHitList()};
    caloHitList.insert(caloHitList.begin(), isolatedHitList.begin(), isolatedHitList.end());

    if (!caloHitList.empty())
    {
        const HitType refView{caloHitList.front()->GetHitType()};
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const HitType view{pCaloHit->GetHitType()};
            const LArCaloHit *const pLArCaloHit{dynamic_cast<const LArCaloHit*>(pCaloHit)};
            if (view != refView || !pLArCaloHit)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
            const unsigned int tpc{pLArCaloHit->GetLArTPCVolumeId()};
            const unsigned int child{pLArCaloHit->GetDaughterVolumeId()};
            if (tpc != m_tpc || child != m_child)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            this->GetLocalCoordinate(pLArCaloHit, localCoords);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TpcHitVolume::Contains(const CaloHit *const pCaloHit) const
{
    const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
    if (!pLArCaloHit)
        return false;

    // ATTN - This assumes only two TPC volumes - ideally, edit in larpandora to ensure consistency
    const unsigned int adjustedDaughterVolumeId{pLArCaloHit->GetLArTPCVolumeId() == 0 ? 2 * pLArCaloHit->GetDaughterVolumeId() :
        1 + 2 * pLArCaloHit->GetDaughterVolumeId()};

    // Might also want to check bounds due to potential cryostat degeneracy
    if (pLArCaloHit->GetLArTPCVolumeId() == m_tpc && adjustedDaughterVolumeId == m_child)
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::GetLocalCoordinate(const CaloHit *pCaloHit, CartesianPointVector &localCoords) const
{
    const HitType view{pCaloHit->GetHitType()};
    const CartesianVector &min{view == TPC_VIEW_U ? m_uMin : view == TPC_VIEW_V ? m_vMin : m_min};
    const CartesianVector &max{view == TPC_VIEW_U ? m_uMax : view == TPC_VIEW_V ? m_vMax : m_max};
    const CartesianVector pos{pCaloHit->GetPositionVector()};
    const double x{pos.GetX()}, z{pos.GetZ()};
    const double xp{(x - min.GetX()) / (max.GetX() - min.GetX())};
    const double zp{(z - min.GetZ()) / (max.GetZ() - min.GetZ())};
    localCoords.emplace_back(CartesianVector(xp, 0., zp));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::InitialiseView(const HitType view)
{
    if (m_viewToCaloHitMap.find(view) == m_viewToCaloHitMap.end())
        m_viewToCaloHitMap[view] = CaloHitList();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::_Echo(const CaloHit *const pCaloHit) const
{
    const HitType view{pCaloHit->GetHitType()};
    const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
    const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
    if (view == TPC_VIEW_U)
    {
        const unsigned int adjustedDaughterVolumeId{pLArCaloHit->GetLArTPCVolumeId() == 0 ? 2 * pLArCaloHit->GetDaughterVolumeId() :
            1 + 2 * pLArCaloHit->GetDaughterVolumeId()};
        std::cout << "U Hit [" << pLArCaloHit->GetLArTPCVolumeId() << "," << pLArCaloHit->GetDaughterVolumeId() << "] => [" <<
            pLArCaloHit->GetLArTPCVolumeId() << "," << adjustedDaughterVolumeId << "]: (" << hitPosition.GetX() <<
            "," << hitPosition.GetZ() << ") TPC " << m_tpc << " Child " << m_child <<
            " Min (" << m_uMin.GetX() << "," << m_uMin.GetZ() <<
            ") Max (" << m_uMax.GetX() << "," << m_uMax.GetZ() << ")" << std::endl;
    }
    if (view == TPC_VIEW_V)
    {
        const unsigned int adjustedDaughterVolumeId{pLArCaloHit->GetLArTPCVolumeId() == 0 ? 2 * pLArCaloHit->GetDaughterVolumeId() :
            1 + 2 * pLArCaloHit->GetDaughterVolumeId()};
        std::cout << "V Hit [" << pLArCaloHit->GetLArTPCVolumeId() << "," << pLArCaloHit->GetDaughterVolumeId() << "] => [" <<
            pLArCaloHit->GetLArTPCVolumeId() << "," << adjustedDaughterVolumeId << "]: (" << hitPosition.GetX() <<
            "," << hitPosition.GetZ() << ") TPC " << m_tpc << " Child " << m_child <<
            " Min (" << m_vMin.GetX() << "," << m_vMin.GetZ() <<
            ") Max (" << m_vMax.GetX() << "," << m_vMax.GetZ() << ")" << std::endl;
    }
    if (view == TPC_VIEW_W)
    {
        const unsigned int adjustedDaughterVolumeId{pLArCaloHit->GetLArTPCVolumeId() == 0 ? 2 * pLArCaloHit->GetDaughterVolumeId() :
            1 + 2 * pLArCaloHit->GetDaughterVolumeId()};
        std::cout << "W Hit [" << pLArCaloHit->GetLArTPCVolumeId() << "," << pLArCaloHit->GetDaughterVolumeId() << "] => [" <<
            pLArCaloHit->GetLArTPCVolumeId() << "," << adjustedDaughterVolumeId << "]: (" << hitPosition.GetX() <<
            "," << hitPosition.GetZ() << ") TPC " << m_tpc << " Child " << m_child <<
            " Min (" << m_min.GetX() << "," << m_min.GetY() << "," << m_min.GetZ() <<
            ") Max (" << m_max.GetX() << "," << m_max.GetY() << "," << m_max.GetZ() << ")" << std::endl;
    }
}

} // namespace lar_content

