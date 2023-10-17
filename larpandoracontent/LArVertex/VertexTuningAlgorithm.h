/**
 *  @file   larpandoracontent/LArVertex/VertexTuningAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_TUNING_ALGORITHM_H
#define LAR_VERTEX_TUNING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexTuningAlgorithm class
 */
class VertexTuningAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexTuningAlgorithm();

private:
    pandora::StatusCode Run();
    void Refine(const pandora::CaloHitList &caloHitList, pandora::CartesianVector &vertex) const;
    float GetScore(const pandora::CaloHitVector &caloHitVector, const pandora::CartesianVector &vertex) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;  ///< The name of the calo hit list
    std::string m_vertexListName;   ///< The name of the vertex list
    std::string m_outputVertexListName;   ///< The name of the output vertex list
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_TUNING_ALGORITHM_H
