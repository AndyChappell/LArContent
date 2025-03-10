/**
 *  @file   larpandoradlcontent/LArVertexing/DlSecondaryVertexingAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SECONDARY_VERTEXING_ALGORITHM_H
#define LAR_DL_SECONDARY_VERTEXING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArVertex/DlVertexingBaseAlgorithm.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningTrackShowerIdAlgorithm class
 */
class DlSecondaryVertexingAlgorithm : public DlVertexingBaseAlgorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DlSecondaryVertexingAlgorithm();

    ~DlSecondaryVertexingAlgorithm();

private:
    class Canvas
    {
    public:
        /**
         *  @brief Default constructor
         */
        Canvas(const pandora::HitType view, const int width, const int height, const int colOffset, const int rowOffset, const float xMin,
            const float xMax, const float zMin, const float zMax);

        virtual ~Canvas();

        pandora::HitType m_view;
        float **m_canvas;
        bool **m_visited;
        const int m_width;
        const int m_height;
        const int m_colOffset;
        const int m_rowOffset;
        const float m_xMin;
        const float m_xMax;
        const float m_zMin;
        const float m_zMax;
    };

    typedef std::map<pandora::HitType, Canvas *> CanvasViewMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /*
     *  @brief  Create input for the network from a calo hit list
     *
     *  @param  caloHits The CaloHitList from which the input should be made
     *  @param  view The wire plane view
     *  @param  xMin The minimum x coordinate for the hits
     *  @param  xMax The maximum x coordinate for the hits
     *  @param  zMin The minimum x coordinate for the hits
     *  @param  zMax The maximum x coordinate for the hits
     *  @param  networkInput The TorchInput object to populate
     *  @param  pixelVector The output vector of populated pixels
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode MakeNetworkInputFromHits(const pandora::CaloHitList &caloHits, const pandora::HitType view, const float xMin,
        const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const;

    /*
     *  @brief  Create a list of vertices from canvases
     *
     *  @param  canvases The input canvases
     *  @param  positionVector The output vector of wire plane positions
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode GetNetworkVertices(const CanvasViewMap &canvases, pandora::CartesianPointVector &positionVector) const;

    /*
     *  @brief  Create a list of vertices from a canvas
     *
     *  @param  canvases The input canvases
     *  @param  positionVector The output vector of wire plane positions
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode GetVerticesFromCanvas(const Canvas &canvas, pandora::CartesianPointVector &vertices) const;

    /**
     *  @brief  Determine if the pixel under consideration is part of a peak and grow that peak to include all connected pixels of equal value
     *
     *  @param  canvas The canvas within which peaks are sought
     *  @param  col The column of the pixel under consideration
     *  @param  row The row of the pixel under consideration
     *  @param  intensity The target intensity of the candidate peak
     *  @param  peak The output vector of pixels constituting the peak under consideration
     *
     *  @return true if we found a better peak while growing the current region, false otherwise
     */
    bool GrowPeak(const Canvas &canvas, int col, int row, float intensity, std::vector<std::pair<int, int>> &peak) const;

    /**
     *  @brief Create a vertex list from the candidate vertices.
     *
     *  @param  candidates The candidate positions with which to create the list.
     *
     *  @return The StatusCode resulting from the function
     */
    pandora::StatusCode MakeCandidateVertexList(const pandora::CartesianPointVector &positions);

    int m_event;                ///< The current event number
    bool m_visualise;           ///< Whether or not to visualise the candidate vertices
    bool m_writeTree;           ///< Whether or not to write validation details to a ROOT tree
    std::string m_rootTreeName; ///< The ROOT tree name
    std::string m_rootFileName; ///< The ROOT file name
    std::mt19937 m_rng;         ///< The random number generator
};

} // namespace lar_dl_content

#endif // LAR_DL_SECONDARY_VERTEXING_ALGORITHM_H
