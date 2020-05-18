//
// Created by Sylvain Renevey on 5/11/18.
//

#ifndef CPP_PROPAGATOR_EOM_H
#define CPP_PROPAGATOR_EOM_H

#include "SatOps/dimensions.h"
#include "SatOps/Vector.h"

class EnvironmentModel;
class Spacecraft;
class StateVector;

/** Class defining the equations of motion. */
class EqMotion {
public:
    /**
    *  @param sc           Instance of the Spacecraft object containing the properties of the spacecraft.
    *  @param env_model    Instance of the EnvironmentModel class containing the parameters of the environment model.
    *  @param initial_et   Initial ephemeris time.
    */
    EqMotion(const Spacecraft &sc, const EnvironmentModel &env_model, double initial_et);

    /** Computes the time derivative of the state at a given time.
     *
     *  This function computes the time derivative of the thirteen states of the spacecraft. It is used by boost's integrator to propagate the states.
     *
     *  @param[in] state                Inertial state of the spacecraft.
     *  @param[in,out] state_derivative     State vector where the time derivatives of the states are stored.
     *  @param[in,out] t                    Time at which the time derivatives are taken
     */
    void operator()(StateVector &state, StateVector &state_derivative, double t);

    /** Computes the gravitational perturbation induced by the nonspherical geopotential.
     *
     * This function computes the perturbations produced by the nonspherical geopotential.
     * The algorithm used to compute the spherical harmonics expansion is described in Montenbruck, O., Gill, E.,
     * <em>Satellites Orbits</em>, Springer-Verlag Berlin Heidelberg, 2000. This method uses the EGM2008 geopotential model.
     *
     *  @param[out] acc              Reference to the acceleration vector.
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state			State vector of the spacecraft.
     *  @param[in] et				Ephemeris time.
     */
    void geopotential(Vector3<dimension::Acceleration> &acc, Vector3d &torque, const StateVector &state, double et) const;

    /** Computes the perturbations induced by the drag.
     *
     *  This function computes the perturbations produced by the atmospheric drag. Based on the settings of the environment model, the atmospheric density
     *  is retrieved either from an exponential model or from NASA's EarthGRAM2016.
     *
     *  @param[out] acc              Reference to the acceleration vector.
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state			State vector of the spacecraft expressed in the J2000 frame (ECI).
     *  @param[in] et				Ephemeris time.
     *  @param[in] elapsed_time     Time elapsed since the beginning of the integration (sec).
     */
    void drag(Vector3<dimension::Acceleration> &acc, Vector3d &torque, const StateVector &state, double et, double elapsed_time) const;

    /** Computes the perturbations induced by the solar radiation pressure.
     *
     *  @param[out] acc              Reference to the acceleration vector.
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state		    State vector of the spacecraft.
     *  @param[in] et				Ephemeris time.
     */
    void solarRadiationPressure(Vector3<dimension::Acceleration> &acc, Vector3d &torque, const StateVector& state, double et) const;

    /** Computes the perturbations induced by the Earth magnetic field.
     *
     *  @param[out] torque           Reference to the torque vector.
     *  @param[in] state			State vector of the spacecraft expressed in the J2000 frame.
     *  @param[in] et				Ephemeris time.
     */
    void magneticPerturbations(Vector3d &torque, const StateVector &state, double et) const;

    /** Computes the perturbations induced by third bodies.
     *
     *  @param[out] acc             Reference to the acceleration vector.
     *  @param[in] state			State vector of the spacecraft.
     *  @param[in] et				Ephemeris time.
     */
    void thirdBodyEffect(Vector3<dimension::Acceleration> &acc, const StateVector &state, double et) const;

private:
    const EnvironmentModel &m_env_model;
    double m_initial_et;
    const Spacecraft &m_spacecraft;
};

#endif //CPP_PROPAGATOR_EOM_H
