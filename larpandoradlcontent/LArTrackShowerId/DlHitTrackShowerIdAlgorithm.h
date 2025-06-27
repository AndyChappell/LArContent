/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlHitTrackShowerIdAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{

/**
 *  @brief  DlHitTrackShowerIdAlgorithm class
 */
class DlHitTrackShowerIdAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlHitTrackShowerIdAlgorithm();

    virtual ~DlHitTrackShowerIdAlgorithm();

private:
    enum Category
    {
        MIP = 1,
        HIP = 2,
        SHOWER = 3,
        MICHEL = 4,
        DIFFUSE = 5
    };

    pandora::StatusCode Run();

    /**
     *  @brief  Prepare training sample for the network
     *
     *  @return The StatusCode indicating success or failure in producing a training sample
     */
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief  Run network inference
     *
     *  @return The StatusCode indicating success or failure in running inference
     */
    pandora::StatusCode Infer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< Name of input calo hit list
    bool m_trainingMode; ///< Training mode
};

} // namespace lar_dl_content

#endif // LAR_DL_HIT_TRACK_SHOWER_ID_ALGORITHM_H
