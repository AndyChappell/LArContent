/**
 *  @file   larpandoracontent/LArUtility/ElectronEnergyAlgorithm.h
 *
 *  @brief  Header file for the electron energy algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ELECTRON_ENERGY_ALGORITHM_H
#define LAR_ELECTRON_ENERGY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ElectronEnergyAlgorithm class
 */
class ElectronEnergyAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PlotGroundTruth() const;
    pandora::StatusCode OutputEnergies(const unsigned int targetPdg = 11) const;
    pandora::StatusCode VisualizeByParent(const unsigned int targetPdg = 11) const;

    std::string m_outputFilename;   ///< The output filename
};

} // namespace lar_content

#endif // #ifndef LAR_ELECTRON_ENERGY_ALGORITHM_H

