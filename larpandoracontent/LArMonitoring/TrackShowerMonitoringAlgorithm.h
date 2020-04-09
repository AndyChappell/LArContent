/**
 *  @file   larpandoracontent/LArMonitoring/TrackShowerMonitoringAlgorithm.h
 *
 *  @brief  Header file for the track shower monitoring algorithm class
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_SHOWER_MONITORING_ALGORITHM_H
#define LAR_TRACK_SHOWER_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief TrackShowerMonitoringAlgorithm class
 */
class TrackShowerMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackShowerMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Visualize PFOs
     *
     *  @param  listName the PFO list name
     */
    void VisualizePfoList(const std::string &listName) const;

    /**
     *  @brief  Visualize available clusters
     *
     *  @param  listName the cluster list name
     */
    void VisualizeAvailableClusterList(const std::string &listName) const;

    /**
     *  @brief  Visualize the calo hit truth
     */
    void VisualizeCaloHitTruth() const;

    bool                    m_showTruth;                ///< Whether to show calo hits with truth info

    pandora::StringVector   m_clusterListNames;         ///< Names of cluster lists to show
    pandora::StringVector   m_pfoListNames;             ///< Names of pfo lists to show
    pandora::StringVector   m_caloHitListNames;         ///< Names of calo hit lists to show
    std::string             m_caloHitList2DName;        ///< Name of 2D calo hit list to show

    std::map<pandora::HitType, std::string> m_viewToNameMap;
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_SHOWER_MONITORING_ALGORITHM_H
