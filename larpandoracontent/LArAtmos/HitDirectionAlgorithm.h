/**
 *  @file   larpandoracontent/LArMonitoring/HitDirectionAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_DIRECTION_ALGORITHM_H
#define LAR_HIT_DIRECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  HitDirectionAlgorithm class
 */
class HitDirectionAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    HitDirectionAlgorithm();

    virtual ~HitDirectionAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_writeFile;               // Whether to produce ROOT output file
    std::string m_filename;         // The filename of the ROOT output file
    std::string m_treename;         // The name of the ROOT tree
};

} // namespace lar_content

#endif // LAR_HIT_DIRECTION_ALGORITHM_H
