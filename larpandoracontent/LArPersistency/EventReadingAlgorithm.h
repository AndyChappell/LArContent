/**
 *  @file   larpandoracontent/LArPersistency/EventReadingAlgorithm.h
 *
 *  @brief  Header file for the event reading algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_READING_ALGORITHM_H
#define LAR_EVENT_READING_ALGORITHM_H 1

#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "Pandora/PandoraInputTypes.h"

#include "Persistency/PandoraIO.h"

namespace pandora
{
class FileReader;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  EventReadingAlgorithm class
 */
class EventReadingAlgorithm : public pandora::ExternallyConfiguredAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventReadingAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventReadingAlgorithm();

    /**
     *  @brief  External event reading parameters class
     */
    class ExternalEventReadingParameters : public pandora::ExternalParameters
    {
    public:
        std::string        m_geometryFileName;  ///< Name of the file containing geometry information
        std::string        m_eventFileNameList; ///< Colon-separated list of file names to be processed
        pandora::InputUInt m_skipToEvent;       ///< Index of first event to consider in input file
    };

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();

    void MoveToNextEventFile();

    pandora::StatusCode ReplaceEventFileReader(const std::string &fileName);

    pandora::FileType GetFileType(const std::string &fileName) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string               m_geometryFileName;      ///< Name of the file containing geometry information
    std::string               m_eventFileName;         ///< Name of the current file containing event information
    pandora::StringVector     m_eventFileNameVector;   ///< Vector of file names to be processed

    unsigned int m_skipToEvent;      ///< Index of first event to consider in first input file
    bool         m_useLArCaloHits;   ///< Whether to read lar calo hits, or standard pandora calo hits
    bool         m_useLArMCParticles;///< Whether to read lar mc particles, or standard pandora mc particles

    pandora::FileReader *m_pEventFileReader; ///< Address of the event file reader
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_READING_ALGORITHM_H
