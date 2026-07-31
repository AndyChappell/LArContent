/**
 *  @file   larpandoracontent/LArPersistency/EventWritingAlgorithm.h
 *
 *  @brief  Header file for the event writing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_WRITING_ALGORITHM_H
#define LAR_EVENT_WRITING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Persistency/PandoraIO.h"

namespace pandora
{
class FileWriter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  EventWritingAlgorithm class
 */
class EventWritingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventWritingAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventWritingAlgorithm();

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();

    bool PassNuanceCodeFilter() const;
    bool PassMCParticleFilter() const;
    bool PassNeutrinoVertexFilter() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::FileType m_geometryFileType;
    pandora::FileType m_eventFileType;

    pandora::FileWriter *m_pEventFileWriter;
    pandora::FileWriter *m_pGeometryFileWriter;

    bool m_shouldWriteGeometry;
    bool m_writtenGeometry;
    std::string m_geometryFileName;

    bool m_shouldWriteEvents;
    std::string m_eventFileName;

    bool m_writtenEventGlobalHeader;

    bool m_shouldWriteMCRelationships;
    bool m_shouldWriteTrackRelationships;

    bool m_shouldOverwriteEventFile;
    bool m_shouldOverwriteGeometryFile;

    bool m_useLArCaloHits;
    bool m_useLArMCParticles;

    bool m_shouldFilterByNuanceCode;
    int  m_filterNuanceCode;

    bool         m_shouldFilterByMCParticles;
    bool         m_neutrinoInducedOnly;
    unsigned int m_matchingMinPrimaryHits;
    unsigned int m_nNonNeutrons;
    unsigned int m_nMuons;
    unsigned int m_nElectrons;
    unsigned int m_nProtons;
    unsigned int m_nPhotons;
    unsigned int m_nChargedPions;

    bool  m_shouldFilterByNeutrinoVertex;
    float m_detectorHalfLengthX;
    float m_detectorHalfLengthY;
    float m_detectorHalfLengthZ;
    float m_coordinateOffsetX;
    float m_coordinateOffsetY;
    float m_coordinateOffsetZ;
    float m_selectedBorderX;
    float m_selectedBorderY;
    float m_selectedBorderZ;
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_WRITING_ALGORITHM_H
