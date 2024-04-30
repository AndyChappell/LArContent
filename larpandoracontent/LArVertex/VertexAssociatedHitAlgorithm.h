/**
 *  @file   larpandoracontent/LArTwoDReco/LArAssociatedHit/VertexAssociatedHitAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_ASSOCIATED_HIT_ALGORITHM_H
#define LAR_VERTEX_ASSOCIATED_HIT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  VertexAssociatedHitAlgorithm class
 */
class VertexAssociatedHitAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexAssociatedHitAlgorithm();

private:
    pandora::StatusCode Run();
    void IdentifyAssociatedHits(const pandora::CaloHitList &caloHitList, const pandora::VertexList &vertexList) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;      ///< The name of the calo hit list
    std::string m_vertexListName;       ///< The name of the vertex list
    std::string m_vetoedHitListName;    ///< The name of the output list of vetoed hits
    std::string m_retainedHitListName;  ///< The name of the output list of retained hits
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_ASSOCIATED_HIT_ALGORITHM_H
