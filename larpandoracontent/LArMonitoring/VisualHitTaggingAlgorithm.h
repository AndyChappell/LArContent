/**
 *  @file   larpandoracontent/LArMonitoring/VisualHitTaggingAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_VISUAL_HIT_TAGGING_ALGORITHM_H
#define LAR_VISUAL_HIT_TAGGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VisualHitTaggingAlgorithm class
 */
class VisualHitTaggingAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VisualHitTaggingAlgorithm();
    virtual ~VisualHitTaggingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_caloHitListNames;    ///< Name of input calo hit lists
    bool m_tagPid;
    bool m_tagPrimary;
    bool m_retainNonReconstructable;
};

} // namespace lar_content

#endif // LAR_VISUAL_HIT_TAGGING_ALGORITHM_H
