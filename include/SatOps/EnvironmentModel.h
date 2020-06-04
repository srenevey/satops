//
// Created by Sylvain Renevey on 1/26/18.
//

#ifndef SATOPS_ENVIRONMENTMODEL_H
#define SATOPS_ENVIRONMENTMODEL_H

#include <array>
#include <cmath>
#include <filesystem>
#include <memory>
#include <vector>
#include "Body.h"
#include "dimensions.h"

#include "ReferenceFrame.h"
#include "Vector.h"

class Atmosphere;
class GravitationalField;
class MagneticField;
class StateVector;

/** Environment model used for the simulation. */
class EnvironmentModel {
public:
    EnvironmentModel();

    /**
     * @param central_body          Central body
     * @param gravity_field         Pointer to a gravitational field.
     * @param magnetic_field        Pointer to a magnetic field.
     * @param atmosphere            Pointer to an atmospheric model.
     * @param third_bodies          List of celestial bodies that are accounted for
     * @param srp_flag              Flag indicating if the solar radiation pressure is accounted for
     */
    explicit EnvironmentModel(
            Body central_body,
            GravitationalField* gravity_field = nullptr,
            MagneticField* magnetic_field = nullptr,
            Atmosphere* atmosphere = nullptr,
            std::vector<Body> third_bodies = {},
            bool srp_flag = false);

    [[nodiscard]] Body centralBody() const;
    [[nodiscard]] std::vector<Body> thirdBodies() const;
    [[nodiscard]] bool isDrag() const;
    [[nodiscard]] bool isSrp() const;
    [[nodiscard]] bool isGravityField() const;
    [[nodiscard]] bool isMagField() const;


    /** Retrieves the density from the atmosphere model.
     *
     *  @param state                State vector of the spacecraft.
     *  @param elapsed_time         Time elapsed since the beginning of the integration (sec)
     *  @return                     Atmosphere density in kg/km<sup>3</sup>.
     */
    [[nodiscard]] double getAtmDensity(const StateVector& state, double elapsed_time) const;


    /** Computes the position vector from the central body to the given body.
     *
     *  @param body         Body to retrieve the position vector for.
     *  @param et			Ephemeris time.
     *  @param frame        Reference frame in which the vector is expressed.
     */
    [[nodiscard]] Vector3<dimension::Distance> getBodyPosition(const Body& body, double et, ReferenceFrame frame) const;

    /** Computes the position vector from the sun to the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Position vector expressed in the reference frame of the state vector.
     */
    [[nodiscard]] Vector3<dimension::Distance> sunSpacecraftVector(const StateVector& state) const;

    /** Computes the shadow condition on the spacecraft.
     *
     *  @param state        State vector of the spacecraft.
     *  @param et			Ephemeris time.
     *  @return             Decimal number with value 0 if no shadow, 1 if in umbra, and between 0 and 1 if in penumbra.
     */
    [[nodiscard]] double isInShadow(const StateVector& state) const;

    /** Computes the acceleration resulting from the gravitational field.
     *
     * @param state     State vector of the spacecraft.
     * @return          Acceleration vector in km/s<sup>2</sup> expressed in the reference frame of the state vector.
     */
    [[nodiscard]] Vector3d computeGravAcc(const StateVector& state) const;


    /** Computes the Earth magnetic field at the given position. The computation is based on a dipole model.
     *
     *  @param state    State vector of the spacecraft.
     *  @param et		Ephemeris time.
     *  @return         Local magnetic field vector in Tesla, expressed in the reference frame of the state vector.
     */
    [[nodiscard]] Vector3d computeMagneticField(const StateVector& state) const;

private:
    Body m_central_body; /*!< Central body. */
    GravitationalField* m_gravity_field; /*!< Pointer to the gravitational field model. */
    MagneticField* m_magnetic_field; /*!< Pointer to the magnetic field model. */
    Atmosphere* m_atmosphere; /*!< Pointer to the atmospheric model. */
    std::vector<Body> m_third_bodies; /*!< Vector containing the celestial bodies accounted for in the computation of third body perturbations. */
    bool m_srp_flag; /*!< Flag indicating if the perturbations due to solar pressure radiation are active. */
};



#endif //SATOPS_ENVIRONMENTMODEL_H
