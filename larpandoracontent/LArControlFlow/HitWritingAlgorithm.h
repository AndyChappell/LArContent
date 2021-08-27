/**
 *  @file   larpandoracontent/LArControlFlow/HitWritingAlgorithm.h
 *
 *  @brief  Header file for the pre processing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_WRITING_ALGORITHM_H
#define LAR_HIT_WRITING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  HitWritingAlgorithm class
 */
class HitWritingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HitWritingAlgorithm();

    ~HitWritingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_HIT_WRITING_ALGORITHM_H
