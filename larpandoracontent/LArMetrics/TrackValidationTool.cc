/**
 *  @file   larpandoracontent/LArMetrics/TrackValidationTool.cc
 *
 *  @brief  Implementation of the track validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMetrics/TrackValidationTool.h"

using namespace pandora;

namespace lar_content
{

TrackValidationTool::TrackValidationTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const /*pMCNu*/, 
    const LArHierarchyHelper::MCMatchesVector &/*mcMatchesVec*/, const MCParticleVector &targetMC, 
    const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    TrackTreeVars trackTreeVars;
    trackTreeVars.m_run = this->GetPandora().GetRun();
    trackTreeVars.m_subrun = this->GetPandora().GetSubrun();
    trackTreeVars.m_event = this->GetPandora().GetEvent();

    this->MichelValidation(pAlgorithm, targetMC, bestRecoMatch, trackTreeVars);
    this->FillTree(trackTreeVars);
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::MichelValidation(const Algorithm *const pAlgorithm, const MCParticleVector &targetMC, const PfoVector &bestRecoMatch,
    TrackTreeVars &trackTreeVars)
{
    for (unsigned int iMC = 0; iMC < targetMC.size(); ++iMC)
    {
        const MCParticle *const pMCParent(targetMC.at(iMC));
        const MCParticle *pMCMichel(nullptr);
        int this_hasMichel(0), this_hasTargetMichel(0), this_hasRecoMichel(0);
        int this_michelIndex(-1), this_michelIsChild(0), this_michelIsShower(0);

        // Find MC michel
        if (std::abs(pMCParent->GetParticleId()) == 13)
        {
            for (const MCParticle *const pMCChild : pMCParent->GetDaughterList())
            {
                if (std::abs(pMCChild->GetParticleId()) != 11)
                    continue;

                if ((pMCParent->GetEndpoint() - pMCChild->GetVertex()).GetMagnitude() < 3.f)
                    continue;

                pMCMichel = pMCChild;
                break;
            }
        }

        if (pMCMichel)
        {
            this_hasMichel = 1;

            // Is Michel a reco target?
            const auto michelIter(std::find(targetMC.begin(), targetMC.end(), pMCMichel));
            this_hasTargetMichel = (michelIter != targetMC.end());

            if (this_hasTargetMichel)
            {
                this_michelIndex = michelIter - targetMC.begin();

                // Has target michel been reco'd
                if (bestRecoMatch.at(this_michelIndex))
                {
                    this_hasRecoMichel = 1;
                    this_michelIsShower = LArPfoHelper::IsShower(bestRecoMatch.at(this_michelIndex));

                    // Is parent-child link correct?
                    if (bestRecoMatch.at(iMC))
                    {
                        const PfoList &muonChildren(bestRecoMatch.at(iMC)->GetDaughterPfoList());
                        this_michelIsChild = (std::find(muonChildren.begin(), muonChildren.end(),
                                                        bestRecoMatch.at(this_michelIndex)) != muonChildren.end());
                    }
                }
            }
        }

        trackTreeVars.m_hasMichel.push_back(this_hasMichel);
        trackTreeVars.m_hasTargetMichel.push_back(this_hasTargetMichel);
        trackTreeVars.m_michelIndex.push_back(this_michelIndex);
        trackTreeVars.m_hasRecoMichel.push_back(this_hasRecoMichel);
        trackTreeVars.m_michelIsChild.push_back(this_michelIsChild);
        trackTreeVars.m_michelIsShower.push_back(this_michelIsShower);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::FillTree(TrackTreeVars &trackTreeVars)
{
    IntVector& hasMichel = trackTreeVars.m_hasMichel;
    IntVector& hasTargetMichel = trackTreeVars.m_hasTargetMichel;
    IntVector& hasRecoMichel = trackTreeVars.m_hasRecoMichel;
    IntVector& michelIndex = trackTreeVars.m_michelIndex;
    IntVector& michelIsChild = trackTreeVars.m_michelIsChild;
    IntVector& michelIsShower = trackTreeVars.m_michelIsShower;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "Run", trackTreeVars.m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "Subrun", trackTreeVars.m_subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "Event", trackTreeVars.m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_HasMichel", &hasMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_HasTargetMichel", &hasTargetMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_IsMichelRecod", &hasRecoMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_MichelIndex", &michelIndex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_MichelIsChild", &michelIsChild));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_MichelIsShower", &michelIsShower));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "TrackTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackValidationTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
