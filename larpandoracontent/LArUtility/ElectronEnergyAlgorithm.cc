/**
 *  @file   larpandoracontent/LArUtility/ElectronEnergyAlgorithm.cc
 *
 *  @brief  Implementation of the electron energy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ElectronEnergyAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

using namespace pandora;

namespace lar_content
{

StatusCode ElectronEnergyAlgorithm::Run()
{
    StatusCode status{STATUS_CODE_SUCCESS};
/*    status = this->OutputEnergies(13);
    if (status == STATUS_CODE_SUCCESS)
        status = this->VisualizeByParent(13);*/
    status = this->PlotGroundTruth();

    return status;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronEnergyAlgorithm::PlotGroundTruth() const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    
    const MCParticleList* pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
                
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_foldBackHierarchy = false;
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);
    
    std::set<unsigned int> pdgCodes{};
    CaloHitList michelHits, lowEnergyHits, mipHits, hipHits, showerHits, otherHits;
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const float mips{pCaloHit->GetMipEquivalentEnergy()};
        if (mips < 0.33)
            continue;
        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            if (targetMCParticleToHitsMap.find(pMCParticle) != targetMCParticleToHitsMap.end())
            {
                const unsigned int pdg = std::abs(pMCParticle->GetParticleId());
                if (pdg == 11)
                {
                    const bool isMichel{LArMCParticleHelper::IsDescendentOf(pMCParticle, 13)};
                    if (isMichel)
                    {
                        michelHits.push_back(pCaloHit);
                        if (pMCParticle->GetEnergy() < 0.05)
                            lowEnergyHits.push_back(pCaloHit);
                        else
                            showerHits.push_back(pCaloHit);
                    }
                    else if (pMCParticle->GetEnergy() < 0.05)
                    {
                        lowEnergyHits.push_back(pCaloHit);
                    }
                    else
                    {
                        showerHits.push_back(pCaloHit);
                    }
                }
                else if(pdg == 22)
                {
                    if (pMCParticle->GetEnergy() < 0.05)
                        lowEnergyHits.push_back(pCaloHit);
                    else
                        showerHits.push_back(pCaloHit);
                }
                else if (pdg == 13 || pdg == 211)
                {
                    mipHits.push_back(pCaloHit);
                }
                else if (pdg == 2212 || pdg == 321 || pdg > 1e9)
                {
                    if (pdg > 1e9)
                        std::cout << "Found nuclei" << std::endl;
                    hipHits.push_back(pCaloHit);
                }
                else
                {
                    otherHits.push_back(pCaloHit);
                    pdgCodes.emplace(pdg);
                    if (pdg == 2112)
                        std::cout << "Neutron hit" << std::endl;
                }                    
            }
        }
        catch(...)
        {
        }
    }

    if (michelHits.size() > 0)
    {
        std::cout << "Michel length: " << michelHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &michelHits, "Michel", GREEN));
    }
    if (lowEnergyHits.size() > 0)
    {
        std::cout << "Low E length: " << lowEnergyHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &lowEnergyHits, "Low E", BLACK));
    }
    if (mipHits.size() > 0)
    {
        std::cout << "MIP length: " << mipHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &mipHits, "MIP", BLUE));
    }
    if (hipHits.size() > 0)
    {
        std::cout << "HIP length: " << hipHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hipHits, "HIP", MAGENTA));
    }
    if (showerHits.size() > 0)
    {
        std::cout << "Shower length: " << showerHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &showerHits, "Shower", RED));
    }
    if (otherHits.size() > 0)
    {
        std::cout << "Other length: " << otherHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &otherHits, "Other", GRAY));
    }
    std::cout << "Other codes";
    for (const unsigned int pdg : pdgCodes)
        std::cout << " " << pdg;
    std::cout << std::endl;
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronEnergyAlgorithm::OutputEnergies(const unsigned int targetPdg) const
{
    const MCParticleList* pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    LArMvaHelper::MvaFeatureVector featureVector;
    LArMvaHelper::MvaFeatureVector ancestorsVector;

    for (const MCParticle* pMCParticle : *pMCParticleList)
    {
        const unsigned int pdg = std::abs(pMCParticle->GetParticleId());
        if (pdg == targetPdg)
        {
            const float energy{pMCParticle->GetEnergy()};
            const bool isMichel{LArMCParticleHelper::IsDescendentOf(pMCParticle, 13)};
            const bool isEMShowerDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 22)};
            const bool isNeutronDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 2112)};
            featureVector.push_back(static_cast<double>(isMichel));
            featureVector.push_back(static_cast<double>(isEMShowerDaughter));
            featureVector.push_back(static_cast<double>(isNeutronDaughter));
            featureVector.push_back(static_cast<double>(energy));
            LArMvaHelper::ProduceTrainingExample(m_outputFilename, true, featureVector);
            
            std::vector<unsigned int> ancestors;
            LArMCParticleHelper::GetAncestors(pMCParticle, ancestors);
            ancestorsVector.push_back(static_cast<double>(energy));
            for (const unsigned int ancestorPdg : ancestors)
                ancestorsVector.push_back(static_cast<double>(ancestorPdg));
            LArMvaHelper::ProduceTrainingExample("ancestors.csv", true, ancestorsVector);
/*            if (isMichel)
                std::cout << "Michel energy: " << energy << std::endl;
//            if (isNeutronDaughter)
//                std::cout << "Neutron energy: " << energy << std::endl;
            if (isEMShowerDaughter)
                std::cout << "EM daughter energy: " << energy << std::endl;

            const bool isNMesonDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 111) ||
                LArMCParticleHelper::IsDescendentOf(pMCParticle, 130)};
            if (isNMesonDaughter)
                std::cout << "Meson 0 energy: " << energy << std::endl;

            const bool isCMesonDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 211) ||
                LArMCParticleHelper::IsDescendentOf(pMCParticle, 321)};
            if (isCMesonDaughter)
                std::cout << "Meson C energy: " << energy << std::endl;*/
        }
        featureVector.clear();
        ancestorsVector.clear();
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronEnergyAlgorithm::VisualizeByParent(const unsigned int targetPdg) const
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    
    const MCParticleList* pMCParticleList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
                
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_foldBackHierarchy = false;
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList2D, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);
    
    CaloHitList michelHits, neutronDaughterHits, meson_zHits, meson_cHits, emHits;    
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            if (targetMCParticleToHitsMap.find(pMCParticle) != targetMCParticleToHitsMap.end())
            {
                const unsigned int pdg = std::abs(pMCParticle->GetParticleId());
                if (pdg != targetPdg)
                    continue;
                const bool isMichel{LArMCParticleHelper::IsDescendentOf(pMCParticle, 13)};
                if (isMichel)
                    michelHits.push_back(pCaloHit);
                    
                const bool isNeutronDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 2112)};
                if (isNeutronDaughter)
                    neutronDaughterHits.push_back(pCaloHit);
                    
                const bool isNMesonDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 111) ||
                    LArMCParticleHelper::IsDescendentOf(pMCParticle, 130)};
                if (isNMesonDaughter)
                    meson_zHits.push_back(pCaloHit);
                    
                const bool isCMesonDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 211) ||
                    LArMCParticleHelper::IsDescendentOf(pMCParticle, 321)};
                if (isCMesonDaughter)
                    meson_cHits.push_back(pCaloHit);
                
                const bool isEMShowerDaughter{LArMCParticleHelper::IsDescendentOf(pMCParticle, 22)};
                if (isEMShowerDaughter)
                    emHits.push_back(pCaloHit);
            }
        }
        catch(...)
        {
        }
    }

    if (michelHits.size() > 0)
    {
        std::cout << "Michel length: " << michelHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &michelHits, "Michel", GREEN));
    }
    if (neutronDaughterHits.size() > 0)
    {
        std::cout << "Neutron daughter length: " << neutronDaughterHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &neutronDaughterHits, "n daughter", GRAY));
    }
    if (meson_zHits.size() > 0)
    {
        std::cout << "N Meson length: " << meson_zHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &meson_zHits, "meson0", MAGENTA));
    }
    if (meson_cHits.size() > 0)
    {
        std::cout << "C Meson length: " << meson_cHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &meson_cHits, "mesonC", CYAN));
    }
    if (emHits.size() > 0)
    {
        std::cout << "EM length: " << emHits.size() << std::endl;
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &emHits, "EM", RED));
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronEnergyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputFilename",
        m_outputFilename));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
