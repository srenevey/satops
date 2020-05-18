//
// Created by Sylvain Renevey on 1/26/18.
//

#ifndef SATOPS_ENVIRONMENTMODEL_H
#define SATOPS_ENVIRONMENTMODEL_H

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include "AtmModel.h"
#include "Atmod1.h"
#include "Body.h"
#include "dimensions.h"
#include "ReferenceFrame.h"
#include "Vector.h"

class StateVector;

/** Environment model used for the simulation. */
class EnvironmentModel {
public:
    EnvironmentModel();

    /**
     * @param central_body          Central body
     * @param gp_degree             Degree of expansion of the geopotential
     * @param geopot_model_path     Path to the geopotential model
     * @param third_bodies          List of celestial bodies that are accounted for
     * @param atm_model             Atmospheric model used for drag computation. Current options are None, Exponential, and EarthGRAM.
     * @param earthgram_path        Path to the EarthGRAM 2016 root directory
     * @param epoch                 Initial epoch with format "YYYY-MM-DD hh:mm:ss"
     * @param srp_flag              Flag indicating if the solar radiation pressure is accounted for
     * @param mag_flag              Flag indicating if the magnetic perturbations are accounted for
     * @param igrf_model_path       Path to the IGRF geomagnetic field model
     */
    EnvironmentModel(
            Body central_body,
            int gp_degree = 0,
            const std::string& geopot_model_path = "",
            const std::vector<Body>& third_bodies = std::vector<Body> {},
            AtmModel atm_model = AtmModel::None,
            const std::string& earthgram_path = "",
            const std::string& epoch = "",
            bool srp_flag = false,
            bool mag_flag = false,
            const std::string& igrf_model_path = "");

    [[nodiscard]] Body centralBody() const;
    [[nodiscard]] std::vector<Body> thirdBodies() const;
    [[nodiscard]] int gpDegree() const;
    [[nodiscard]] bool isDrag() const;
    [[nodiscard]] bool isSrp() const;
    [[nodiscard]] bool isMagField() const;


    /** Returns the normalized C coefficient for the given degree and order. */
    [[nodiscard]] double cCoeff(int degree, int order) const;
    /** Returns the normalized S coefficient for the given degree and order. */
    [[nodiscard]] double sCoeff(int degree, int order) const;


    /** Retrieves the density from the atmosphere model specified in the EnvironmentModel.
     *
     *  @param state                State vector of the spacecraft.
     *  @param elapsed_time         Time elapsed since the beginning of the integration (sec)
     *  @param et                   Ephemeris time.
     *  @return                     Atmosphere density in kg/km<sup>3</sup>.
     */
    [[nodiscard]] double atmDensity(const StateVector& state, double elapsed_time, double et) const;

    /** Computes the V and W coefficients required to compute the acceleration.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et           Ephemeris time.
     */
    [[nodiscard]] std::vector<std::array<double, 2>> geopotHarmonicCoeff(const StateVector& state, double et) const;


    /** Computes the position vector from the central body to the given body.
     *
     *  @param body         Body to retrieve the position vector for.
     *  @param et			Ephemeris time.
     *  @param frame        Reference frame in which the vector is expressed.
     */
    [[nodiscard]] Vector3<dimension::Distance> bodyVector(const Body& body, double et, ReferenceFrame frame) const;

    /** Computes the position vector from the sun to the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Position vector expressed in the reference frame of the state vector.
     */
    [[nodiscard]] Vector3<dimension::Distance> sunSpacecraftVector(const StateVector& state, double et) const;

    /** Computes the shadow condition on the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Decimal number with value 0 if no shadow, 1 if in umbra, and between 0 and 1 if in penumbra.
     */
    [[nodiscard]] double isInShadow(const StateVector& state, double et) const;


    /** Computes the Earth magnetic field at the given position. The computation is based on a dipole model.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Local magnetic field vector in Tesla, expressed in the reference frame of the state vector.
     */
    [[nodiscard]] Vector3d magneticField(const StateVector& state, double et) const;

private:
    Body m_central_body; /*!< Central body. */
    int m_gp_degree; /*!< Degree of expansion of the Earth geopotential. */
    std::unique_ptr<double[]> m_c_geopot_coeffs; /*!< Array containing the C normalized coefficients used to compute the geopotential effect. */
    std::unique_ptr<double[]> m_s_geopot_coeffs; /*!< Array containing the S normalized coefficients used to compute the geopotential effect. */
    std::vector<Body> m_third_bodies; /*!< Vector containing the celestial bodies accounted for in the computation of third body perturbations. */
    AtmModel m_atm_model; /*!< Variant of the atmospheric model to be used. Current options are None, Exponential, and EarthGRAM. */
    std::unique_ptr<Atm1> m_earth_gram_atm_model; /*!< Pointer to an instance of EarthGRAM 2016 model's Atm1 object. */
    double m_exp_atm_model[28][4]; /*!< Array containing the exponential atmospheric model */
    bool m_srp_flag; /*!< Flag indicating if the perturbations due to solar pressure radiation are active. */
    bool m_mag_flag; /*!< Flag indicating if the magnetic perturbations are active. */
    std::unique_ptr<double[]> m_g_mag_coeffs; /*!< Array containing the g normalized coefficients used to compute the Earth magnetic field. */
    std::unique_ptr<double[]> m_h_mag_coeffs; /*!< Array containing the h normalized coefficients used to compute the Earth magnetic field. */

    // Methods
    void loadGeopotCoeff(int degree, std::string model_file_name);
    void loadGeomagCoeff(std::string model_file_name);
    void initExpModel();
    int getIndex(int degree, int order) const;
};



#endif //SATOPS_ENVIRONMENTMODEL_H
