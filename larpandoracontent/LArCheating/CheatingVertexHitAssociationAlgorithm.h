/**
 *  @file   larpandoracontent/LArCheating/CheatingVertexHitAssociationAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_VERTEX_HIT_ASSOCIATION_ALGORITHM_H
#define LAR_CHEATING_VERTEX_HIT_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingVertexHitAssociationAlgorithm::Algorithm class
 */
class CheatingVertexHitAssociationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingVertexHitAssociationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitListName; ///< The name of the input calo hit list
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_VERTEX_HIT_ASSOCIATION_ALGORITHM_H
