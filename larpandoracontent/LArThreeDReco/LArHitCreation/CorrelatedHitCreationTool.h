/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/CorrelatedHitCreationTool.h
 *
 *  @brief  Header file for the correlated hit creation tool.
 *
 *  $Log: $
 */
#ifndef LAR_CORRELATED_HIT_CREATION_TOOL_H
#define LAR_CORRELATED_HIT_CREATION_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar_content
{

/**
 *  @brief  CorrelatedHitCreationTool class
 */
class CorrelatedHitCreationTool : public HitCreationBaseTool
{
public:
    typedef HitCreationBaseTool::ProtoHit ProtoHit;
    typedef HitCreationBaseTool::ProtoHitVector ProtoHitVector;

    /**
     *  @brief  Default constructor
     */
    CorrelatedHitCreationTool();

    /**
     *  @brief  Destructor
     */
    virtual ~CorrelatedHitCreationTool();

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo the address of the pfo
     *  @param  inputTwoDHits the vector of input two dimensional hits
     *  @param  protoHitVector to receive the new three dimensional proto hits
     */
    void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector);

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CORRELATED_HIT_CREATION_TOOL_H
