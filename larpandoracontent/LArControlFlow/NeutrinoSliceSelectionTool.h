/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoSliceSelectionTool.h
 *
 *  @brief  Header file for the neutrino slice selection tool class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_SLICE_SELECTION_TOOL_H
#define LAR_NEUTRINO_SLICE_SELECTION_TOOL_H 1

#include "larpandoracontent/LArControlFlow/GenericSliceSelectionTool.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoSliceSelectionTool class
 */
class NeutrinoSliceSelectionTool : public GenericSliceSelectionTool
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoSliceSelectionTool();

protected:
    /**
     *  @brief  Template method to determine if an MC particle matches the target criteria for slice selection. Return true if match.
     *
     *  @param  mcParticle the MC particle to check
     */
    bool IsTarget(const pandora::MCParticle *const mcParticle) const;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_SLICE_SELECTION_TOOL_H
