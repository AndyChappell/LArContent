/**
 *  @file   larpandoracontent/LArUtility/TpcHitVolume.cc
 *
 *  @brief  Implementation of the list merging algorithm class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArUtility/TpcHitVolume.h"

using namespace pandora;

namespace lar_content
{

TpcHitVolume::TpcHitVolume(const unsigned int cryostat, const unsigned int tpc, const unsigned int child, const pandora::CartesianVector &center,
    const pandora::CartesianVector &length, const LArTransformationPlugin *const pTransform) : m_cryostat{cryostat}, m_tpc{tpc}, m_child{child},
    m_min{center - length * 0.5}, m_max{center + length * 0.5}, m_uMin{0., 0., 0.}, m_uMax{0., 0., 0.}, m_vMin{0., 0., 0.}, m_vMax{0., 0., 0.}
{
    std::cout << "W(min): " << m_min << std::endl;
    std::cout << "W(max): " << m_max << std::endl;

    m_uMin = CartesianVector(m_min.GetX(), m_min.GetY(), static_cast<float>(pTransform->YZtoU(m_min.GetY(), m_min.GetZ())));
    m_uMax = CartesianVector(m_max.GetX(), m_max.GetY(), static_cast<float>(pTransform->YZtoU(m_max.GetY(), m_max.GetZ())));
    m_vMin = CartesianVector(m_min.GetX(), m_min.GetY(), static_cast<float>(pTransform->YZtoV(m_min.GetY(), m_min.GetZ())));
    m_vMax = CartesianVector(m_max.GetX(), m_max.GetY(), static_cast<float>(pTransform->YZtoV(m_max.GetY(), m_max.GetZ())));

    std::cout << "V(min): " << m_vMin << std::endl;
    std::cout << "V(max): " << m_vMax << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TpcHitVolume::~TpcHitVolume()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::Add(const pandora::CaloHitList &caloHitList)
{
    for (const CaloHit *const pCaloHit : caloHitList)
        this->Add(pCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TpcHitVolume::Add(const pandora::CaloHit *const pCaloHit)
{
    if (this->Contains(pCaloHit))
    {
        const HitType view{pCaloHit->GetHitType()};
        this->InitialiseView(view);
        m_viewToCaloHitMap[view].emplace_back(pCaloHit);

        const CartesianVector &hitPosition{pCaloHit->GetPositionVector()};
        const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
/*        if (view == TPC_VIEW_U)
        {
            const unsigned int adjustedDaughterVolumeId{pLArCaloHit->GetLArTPCVolumeId() == 0 ? 2 * pLArCaloHit->GetDaughterVolumeId() :
                1 + 2 * pLArCaloHit->GetDaughterVolumeId()};
            std::cout << "U Hit [" << pLArCaloHit->GetLArTPCVolumeId() << "," << pLArCaloHit->GetDaughterVolumeId() << "] => [" <<
                pLArCaloHit->GetLArTPCVolumeId() << "," << adjustedDaughterVolumeId << "]: (" << hitPosition.GetX() <<
                "," << hitPosition.GetZ() << ") TPC " << m_tpc << " Child " << m_child <<
                " Min (" << m_uMin.GetX() << "," << m_uMin.GetZ() <<
                ") Max (" << m_uMax.GetX() << "," << m_uMax.GetZ() << ")" << std::endl;
        }*/
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
/*        if (view == TPC_VIEW_W)
        {
            const unsigned int adjustedDaughterVolumeId{pLArCaloHit->GetLArTPCVolumeId() == 0 ? 2 * pLArCaloHit->GetDaughterVolumeId() :
                1 + 2 * pLArCaloHit->GetDaughterVolumeId()};
            std::cout << "W Hit [" << pLArCaloHit->GetLArTPCVolumeId() << "," << pLArCaloHit->GetDaughterVolumeId() << "] => [" <<
                pLArCaloHit->GetLArTPCVolumeId() << "," << adjustedDaughterVolumeId << "]: (" << hitPosition.GetX() <<
                "," << hitPosition.GetZ() << ") TPC " << m_tpc << " Child " << m_child <<
                " Min (" << m_min.GetX() << "," << m_min.GetY() << "," << m_min.GetZ() <<
                ") Max (" << m_max.GetX() << "," << m_max.GetY() << "," << m_max.GetZ() << ")" << std::endl;
        }*/
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TpcHitVolume::Contains(const pandora::CaloHit *const pCaloHit) const
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

void TpcHitVolume::InitialiseView(const pandora::HitType view)
{
    if (m_viewToCaloHitMap.find(view) == m_viewToCaloHitMap.end())
        m_viewToCaloHitMap[view] = CaloHitList();
}

} // namespace lar_content

