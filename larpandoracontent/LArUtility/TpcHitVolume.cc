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
    const pandora::CartesianVector &length) : m_cryostat{cryostat}, m_tpc{tpc}, m_child{child}, m_min{center - length * 0.5f},
    m_max{center + length * 0.5f}
{
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
        std::cout << "Hit: (" << hitPosition.GetX() << "," << hitPosition.GetZ() << ") TPC " << m_tpc << " Child " << m_child <<
            " Min (" << m_min.GetX() << "," << m_min.GetY() << "," << m_min.GetZ() <<
            ") Max (" << m_max.GetX() << "," << m_max.GetY() << "," << m_max.GetZ() << ")" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TpcHitVolume::Contains(const pandora::CaloHit *const pCaloHit) const
{
    const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
    if (!pLArCaloHit)
        return false;

    // Might also want to check bounds due to potential cryostat degeneracy
    if (pLArCaloHit->GetLArTPCVolumeId() == m_tpc && pLArCaloHit->GetDaughterVolumeId() == m_child)
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

