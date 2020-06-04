//
// Created by Sylvain  on 6/2/20.
//

#ifndef SATOPS_ATMOSPHERE_H
#define SATOPS_ATMOSPHERE_H

#include <cmath>
#include <filesystem>
#include <memory>
#include "Atmod1.h"

class StateVector;

/** Atmospheric model of the central body. */
class Atmosphere {
public:
    /** Atmospheric model */
    enum Model {
        earthGRAM, /*!< Earth Global Reference Atmospheric Model 2016. */
        exp /*!< Exponentially decaying atmospheric model. */
    };

    /** Creates an atmospheric model.
     *
     * @param[in] model     Name of the model to use.
     * @param[in] file      Path containing the coefficients of the model.
     * @param[in] epoch     Initial epoch of the simulation (required by EarthGRAM only).
     */
    Atmosphere(Model model, const std::filesystem::path& file, std::string epoch = "");

    /** Retrieves the atmospheric density at the given position.
     *
     * @param[in] state         Current state of the spacecraft.
     * @param[in] elapsed_time  Elapsed time since the start of the simulation.
     * @return Atmospheric density in kg/km<sup>3</sup>.
     */
    double getDensity(const StateVector& state, double elapsed_time) const;

private:
    /** Binary search to find the closest altitude for the exponential model. */
    int searchAltitude(double alt) const;
    void loadExpModel(const std::filesystem::path& file);

    Model m_model;
    std::unique_ptr<Atm1> m_earth_gram_model;
    std::unique_ptr<double[]> m_exp_model_base_alt;
    std::unique_ptr<double[]> m_exp_model_density;
    std::unique_ptr<double[]> m_exp_model_scale_height;
};

#endif //SATOPS_ATMOSPHERE_H
