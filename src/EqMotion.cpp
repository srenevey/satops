//
// Created by Sylvain Renevey on 5/11/18.
//

#include "EqMotion.h"
#include <cmath>
#include <cstddef>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include "transformations.h"
#include "SatOps/EnvironmentModel.h"
#include "SatOps/Quaternion.h"
#include "SatOps/Spacecraft.h"
#include "SatOps/StateVector.h"
extern "C" {
#include "SpiceUsr.h"
};

EqMotion::EqMotion(const Spacecraft &sc, const EnvironmentModel &env_model, double initial_et): m_env_model(env_model), m_initial_et(initial_et), m_spacecraft(sc) {};

void EqMotion::operator()(StateVector &state, StateVector &state_derivative, double t) {

    // Compute the current ephemeris time
    double et = m_initial_et + t;
    state.setTime(et);
    state.normalizeQuaternion();

    // Two-body acceleration
    double factor = -m_env_model.centralBody().mu() / pow(state.position().norm(), 3);
    Vector3<dimension::Acceleration> a_kepler(state.frame(), {factor * state[0], factor * state[1], factor * state[2]});

    // Geopotential
    Vector3<dimension::Acceleration> a_geopot(state.frame());
    Vector3d t_gg(ReferenceFrame::BODY);
    if (m_env_model.gpDegree() > 1) {
        geopotential(a_geopot, t_gg, state, et);
    }

    // Third-body effect
    Vector3<dimension::Acceleration> a_third_body(state.frame());
    if (!m_env_model.thirdBodies().empty()) {
        thirdBodyEffect(a_third_body, state, et);
    }

    // Aerodynamic perturbations
    Vector3<dimension::Acceleration> a_aero(state.frame());
    Vector3d t_aero(ReferenceFrame::BODY);
    if (m_env_model.isDrag() && m_spacecraft.mass() != 0.0) {
        drag(a_aero, t_aero, state, et, t);
    }

    // Solar radiation pressure
    Vector3<dimension::Acceleration> a_srp(state.frame());
    Vector3d t_srp(ReferenceFrame::BODY);
    if (m_env_model.isSrp()) {
        solarRadiationPressure(a_srp, t_srp, state, et);
    }

    // Magnetic perturbations
    Vector3d t_mag(ReferenceFrame::BODY);
    if (m_env_model.isMagField()) {
        magneticPerturbations(t_mag, state, et);
    }

    Vector3<dimension::Acceleration> a_total = a_kepler + a_geopot + a_third_body + a_aero + a_srp;
    Vector3d torques = t_gg + t_aero + t_srp + t_mag;

    state_derivative.setPositionDerivative(state.velocity());
    state_derivative.setVelocityDerivative(a_total);


    // Attitude kinematics
    Quaternion q = state.orientation();
    q.normalize();
    Quaternion dq(0.5 * q.xi() * state.angVelocity());
    state_derivative.setOrientation(dq);

    // Attitude dynamics
    Vector3<dimension::AngularAcceleration> ang_acceleration(m_spacecraft.inertiaMatrix().inverse() * (torques -
                                                                                                       state.angVelocity().cross(
            m_spacecraft.inertiaMatrix() * state.angVelocity())));
    state_derivative.setAngVelocityDerivative(ang_acceleration);
}

void EqMotion::geopotential(Vector3<dimension::Acceleration> &acc, Vector3d &torque, const StateVector &state, double et) const {

    // Creates lambdas to retrieve C and S coefficients from the environment model
    auto C = [&](int degree, int order) { return m_env_model.cCoeff(degree, order); };
    auto S = [&](int degree, int order) { return m_env_model.sCoeff(degree, order); };

    // Computes and retrieve the V and W coefficients
    auto vw_coeffs = m_env_model.geopotHarmonicCoeff(state, et);
    // Computes the getIndex for degree n and order m and creates lambda functions
    auto index = [](int degree, int order) { return degree * (degree + 1) / 2 + order; };
    auto V = [&](int degree, int order) { return vw_coeffs[index(degree, order)][0]; };
    auto W = [&](int degree, int order) { return vw_coeffs[index(degree, order)][1]; };

    // Computes the acceleration
    Vector3<dimension::Acceleration> a_geopot(ReferenceFrame::ITRF93);

    // n = degree, m = order
    for (int n = 0; n < m_env_model.gpDegree() + 1; ++n) {
        for (int m = 0; m < n + 1; ++m) {
            double scale = 0;
            if (m == 0) {
                scale = sqrt((2 * n + 1));
                a_geopot[0] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * (-C(n, 0) * V(n+1, 1));
                a_geopot[1] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * (-C(n, 0) * W(n+1, 1));
            } else {
                scale = sqrt(2 * (2 * n + 1) * boost::math::factorial<double>(n-m) / boost::math::factorial<double>(n+m));
                a_geopot[0] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * 0.5 *
                               ((-C(n, m) * V(n+1, m+1) - S(n, m) * W(n+1, m+1)) +
                                boost::math::factorial<double>(n-m+2) / boost::math::factorial<double>(n-m) *
                                (C(n, m) * V(n+1, m-1) + S(n, m) * W(n+1, m-1)));
                a_geopot[1] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) * 0.5 *
                               ((-C(n, m) * W(n+1, m+1) + S(n, m) * V(n+1, m+1)) +
                                boost::math::factorial<double>(n-m+2) / boost::math::factorial<double>(n-m) *
                                (-C(n, m) * W(n+1, m-1) + S(n, m) * V(n+1, m-1)));
            }
            a_geopot[2] += scale * constants::MU_EARTH_EGM08 / (constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08) *
                           ((n-m+1) * (-C(n, m) * V(n+1, m) - S(n, m) * W(n+1, m)));
        }
    }

    // Rotates back to the original frame. Note: a_geopot already accounts for Coriolis, centripetal, etc. accelerations.
    acc = transformations::rotateToFrame(a_geopot, state.frame(), et);

    // Nadir-pointing vector in the body frame
    Vector3<dimension::Distance> nadir = state.position() / state.position().norm();
    Vector3<dimension::Distance> nadir_bff = transformations::rotateToFrame(nadir, ReferenceFrame::BODY, et,
                                                                            state.orientation());
    torque = 3.0 * constants::MU_EARTH / pow(state.position().norm(), 3) * nadir_bff.cross(
            m_spacecraft.inertiaMatrix() * nadir_bff);
}

void EqMotion::drag(Vector3<dimension::Acceleration> &acc, Vector3d &torque, const StateVector &state, double et, double elapsed_time) const {

    // Retrieves atmospheric density
    double density = m_env_model.atmDensity(state, elapsed_time, et);

    // Computes the velocity relative to Earth's atmosphere
    Vector3<dimension::AngularVelocity> earth_ang_vel(ReferenceFrame::J2000, {0., 0., constants::EARTH_ANGULAR_VELOCITY});
    Vector3<dimension::Velocity> earth_velocity(earth_ang_vel.cross(state.position()));
    Vector3<dimension::Velocity> v_rel = state.velocity() - earth_velocity;
    Vector3<dimension::Velocity> v_rel_bff = transformations::rotateToFrame(v_rel, ReferenceFrame::BODY, et,
                                                                            state.orientation());

    // Force and torque
    Vector3<dimension::Acceleration> acc_bff(ReferenceFrame::BODY);
    for (std::size_t i = 0; i < m_spacecraft.faceNormals().size(); ++i) {
        double cos_theta = m_spacecraft.faceNormals()[i].dot(v_rel_bff) / (v_rel.norm() * v_rel.norm());
        Vector3d force_bff(-0.5 * density * m_spacecraft.dragCoeff() * v_rel.norm() * v_rel_bff *
                                   m_spacecraft.faceAreas()[i] * std::max(cos_theta, 0.));
        acc_bff += force_bff/m_spacecraft.mass();
        torque += m_spacecraft.faceCopPositions()[i].cross(force_bff);
    }
    acc = transformations::rotateToFrame(acc_bff, state.frame(), et, state.orientation());
}

void EqMotion::solarRadiationPressure(Vector3<dimension::Acceleration>& acc, Vector3d& torque, const StateVector& state, double et) const {
    Vector3<dimension::Distance> r_sun2sc = m_env_model.sunSpacecraftVector(state, et);
    double nu = m_env_model.isInShadow(state, et);
    double p_sr = (constants::SOLAR_PRESSURE * 1E-3) * constants::AU_TO_KM / r_sun2sc.norm(); // kg / (km s^2)

    Vector3<dimension::Distance> r_sc2sun_bff = transformations::rotateToFrame((-r_sun2sc) / r_sun2sc.norm(),
                                                                               ReferenceFrame::BODY, et,
                                                                               state.orientation());
    Vector3<dimension::Acceleration> acc_bff(ReferenceFrame::BODY);
    for (std::size_t i = 0; i < m_spacecraft.faceNormals().size(); ++i) {
        double cos_theta = m_spacecraft.faceNormals()[i].dot(r_sc2sun_bff);
        Vector3d force_bff(-nu * p_sr * m_spacecraft.faceAreas()[i] * (2.0 * (m_spacecraft.diffuseReflectionCoeff()[i] / 3.0 + m_spacecraft.specularReflectionCoeff()[i] * cos_theta) * m_spacecraft.faceNormals()[i] + (1.0 - m_spacecraft.specularReflectionCoeff()[i]) * r_sc2sun_bff) * std::max(cos_theta, 0.));
        acc_bff += force_bff/m_spacecraft.mass();
        torque += m_spacecraft.faceCopPositions()[i].cross(force_bff);
    }
    acc = transformations::rotateToFrame(acc_bff, state.frame(), et);
}

void EqMotion::magneticPerturbations(Vector3d& torque, const StateVector& state, double et) const {

    Vector3d B = m_env_model.magneticField(state, et);
    Vector3d B_bff = transformations::rotateToFrame(B, ReferenceFrame::BODY, et, state.orientation());

    torque = m_spacecraft.residualDipole().cross(B_bff);
}

void EqMotion::thirdBodyEffect(Vector3<dimension::Acceleration>& acc, const StateVector& state, double et) const {
    auto bodies = m_env_model.thirdBodies();

    for (auto body: bodies) {
        Vector3<dimension::Distance> r_cb2tb = m_env_model.bodyVector(body, et, state.frame()); // central body to third body
        Vector3<dimension::Distance> r_cb2sc = state.position(); // central body to spacecraft
        Vector3<dimension::Distance> r_sc2tb = r_cb2tb - r_cb2sc; // spacecraft to third body

        double numerator = r_cb2sc.norm() * r_cb2sc.norm() + r_cb2tb.norm() * r_cb2tb.norm() -
                           ((r_cb2tb - r_cb2sc).norm()) * ((r_cb2tb - r_cb2sc).norm());
        double denominator = 2 * r_cb2sc.norm() * r_cb2tb.norm();
        double cos_zeta = numerator / denominator;
        double h = r_cb2sc.norm() / r_cb2tb.norm();
        double B = 0.0;
        double err = 1.0;
        int i = 1;

        while (std::abs(err) > 1E-8 && i < 100) {
            double B_new = B + boost::math::legendre_p(i, cos_zeta) * pow(h, i);
            err = B_new - B;
            B = B_new;
            ++i;
        }

        double beta = 3.0 * B + 3.0 * B * B + B * B * B;
        acc += -body.mu()/pow(r_cb2tb.norm(), 3) * (r_cb2sc - beta * (r_cb2tb - r_cb2sc));
    }
}