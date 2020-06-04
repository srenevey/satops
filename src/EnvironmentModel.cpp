//
// Created by Sylvain Renevey on 1/26/18.
//



#include "SatOps/EnvironmentModel.h"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <utility>
#include "SatOps/Atmosphere.h"
#include "SatOps/constants.h"
#include "SatOps/GravitationalField.h"
#include "SatOps/MagneticField.h"
#include "SatOps/StateVector.h"
#include "transformations.h"
extern "C" {
#include "SpiceUsr.h"
};

EnvironmentModel::EnvironmentModel():
        m_central_body(Body::Earth),
        m_gravity_field(nullptr),
        m_magnetic_field(nullptr),
        m_atmosphere(nullptr),
        m_third_bodies({}),
        m_srp_flag(false)
{}

EnvironmentModel::EnvironmentModel(Body central_body, GravitationalField* gravity_field, MagneticField* magnetic_field, Atmosphere* atmosphere, std::vector<Body> third_bodies, bool srp_flag):
        m_central_body(central_body),
        m_gravity_field(gravity_field),
        m_magnetic_field(magnetic_field),
        m_atmosphere(atmosphere),
        m_third_bodies(std::move(third_bodies)),
        m_srp_flag(srp_flag)
{}

Body EnvironmentModel::centralBody() const {
    return m_central_body;
}

std::vector<Body> EnvironmentModel::thirdBodies() const {
    return m_third_bodies;
}

bool EnvironmentModel::isDrag() const {
    return (m_atmosphere != nullptr);
}

bool EnvironmentModel::isSrp() const {
    return m_srp_flag;
}

bool EnvironmentModel::isGravityField() const {
    return (m_gravity_field != nullptr);
}

bool EnvironmentModel::isMagField() const {
    return (m_magnetic_field != nullptr);
}

double EnvironmentModel::getAtmDensity(const StateVector& state, double elapsed_time) const {
    return m_atmosphere->getDensity(state, elapsed_time);
}

Vector3<dimension::Distance> EnvironmentModel::getBodyPosition(const Body& body, double et, ReferenceFrame frame) const {
    std::string central_body = std::to_string(m_central_body.naifId());
    std::string third_body = std::to_string(body.naifId());

    SpiceDouble third_body_state[6];
    SpiceDouble lt = 0;
    spkezr_c(third_body.c_str(), et, std::string(frame).c_str(), "NONE", central_body.c_str(), third_body_state, &lt);
    return Vector3<dimension::Distance>(frame, third_body_state);
}

Vector3<dimension::Distance> EnvironmentModel::sunSpacecraftVector(const StateVector& state) const {
    Vector3<dimension::Distance> r_sun2sc(state.frame());

    if (m_central_body == Body::Sun) {
        r_sun2sc = state.position();
    }
    else {
        std::string cb = std::to_string(m_central_body.naifId());
        Vector3<dimension::Distance> r_cb2sun = getBodyPosition(Body::Sun, state.time(), state.frame());
        r_sun2sc = -r_cb2sun + state.position();
    }
    return r_sun2sc;
}

double EnvironmentModel::isInShadow(const StateVector& state) const {
    Vector3<dimension::Distance> r_sun2sc = sunSpacecraftVector(state);

    double a = std::asin(constants::R_SUN / r_sun2sc.norm());
    double b = std::asin(m_central_body.radius() / state.position().norm());
    double c = std::acos( state.position().dot(r_sun2sc) / (state.position().norm() * r_sun2sc.norm()));

    double nu = 0.0;
    if (std::abs(a-b) < c && c < a + b) {
        double x = (c * c + a * a - b * b) / (2.0 * c);
        double y = std::sqrt(a * a - x * x);
        double A = a * a * std::acos(x / a) + b * b * std::acos((c - x) / b) - c * y;
        nu = 1.0 - A / (constants::PI * a * a);
    } else if (c >= a + b) {
        nu = 0.0;
    } else if (c < b - a) {
        nu = 1.0;
    }
    return nu;
}

Vector3d EnvironmentModel::computeGravAcc(const StateVector& state) const {
    // Transform the position to ITRF93 (rotating).
    StateVector state_itrf93 = state.transformToFrame(ReferenceFrame::ITRF93);
    Vector3<dimension::Distance> r_itrf93 = state_itrf93.position();
    Vector3d acc = m_gravity_field->getAcceleration(r_itrf93);

    // Transform back to initial frame and return
    return transformations::rotateToFrame(acc, state.frame(), state.time());
}

Vector3d EnvironmentModel::computeMagneticField(const StateVector& state) const {
    //Vector3d m(std::pow(constants::IGRF_2020_RADIUS, 3) * Vector3d(ReferenceFrame::ITRF93, {constants::IGRF_2020_G11, constants::IGRF_2020_H11, constants::IGRF_2020_G01}));

    // Transform the position to ITRF93 (rotating).
    StateVector state_itrf93 = state.transformToFrame(ReferenceFrame::ITRF93);
    Vector3<dimension::Distance> r_itrf93 = state_itrf93.position();
    Vector3d B = m_magnetic_field->getMagneticField(r_itrf93);

    // Compute the magnetic field
    //Vector3d B = 3.0 * (m.dot(r_itrf93) * r_itrf93 - r_itrf93.norm()*r_itrf93.norm()*m) / pow(r_itrf93.norm(), 5);

    // Transform back to initial frame and return
    return transformations::rotateToFrame(B, state.frame(), state.time());
}