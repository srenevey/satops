//
// Created by Sylvain  on 4/12/20.
//

#include "SatOps/StateVector.h"
#include <cstddef>
#include <exception>
#include <iostream>
extern "C" {
#include "SpiceUsr.h"
};

StateVector::StateVector(): m_frame(ReferenceFrame::NONE), m_ephemeris_time(0.0), m_position(Vector3<dimension::Distance>()), m_velocity(Vector3<dimension::Velocity>()), m_orientation(Quaternion()), m_ang_velocity(Vector3<dimension::AngularVelocity>()) {}

StateVector::StateVector(double ephemeris_time, Vector3<dimension::Distance> position, Vector3<dimension::Velocity> velocity, Quaternion orientation, Vector3<dimension::AngularVelocity> ang_velocity) :
        m_ephemeris_time(ephemeris_time),
        m_orientation(orientation)
{
    try {
        if (position.frame() != velocity.frame())
            throw std::invalid_argument("StateVector: The position and velocity must be in the same reference frame.");
        if (ang_velocity.frame() != ReferenceFrame::BODY)
            throw std::invalid_argument("StateVector: The angular velocity must be specified in the body-fixed frame.");

        m_frame = position.frame();
        m_position = position;
        m_velocity = velocity;
        m_ang_velocity = ang_velocity;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

StateVector::StateVector(const std::string& initial_epoch, Vector3<dimension::Distance> position, Vector3<dimension::Velocity> velocity, Quaternion orientation, Vector3<dimension::AngularVelocity> ang_velocity) : StateVector(0., position, velocity, orientation, ang_velocity) {
    str2et_c(initial_epoch.c_str(), &m_ephemeris_time);
}

ReferenceFrame StateVector::frame() const {
    return m_frame;
}

double StateVector::time() const {
    return m_ephemeris_time;
}

Vector3<dimension::Distance> StateVector::position() const {
    return m_position;
}

Vector3<dimension::Velocity> StateVector::velocity() const {
    return m_velocity;
}

Quaternion StateVector::orientation() const {
    return m_orientation;
}

Vector3<dimension::AngularVelocity> StateVector::angVelocity() const {
    return m_ang_velocity;
}

void StateVector::setFrame(ReferenceFrame frame) {
    m_frame = frame;
}

void StateVector::setTime(double time) {
    m_ephemeris_time = time;
}

void StateVector::setPosition(Vector3<dimension::Distance> position) {
    try {
        if (m_frame != position.frame())
            throw std::invalid_argument("StateVector: Cannot set position in the state vector: the reference frames are different.");

        m_position = position;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

void StateVector::setVelocity(Vector3<dimension::Velocity> velocity) {
    try {
        if (m_frame != velocity.frame())
            throw std::invalid_argument("StateVector: Cannot set velocity in the state vector: the reference frames are different.");

        m_velocity = velocity;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

void StateVector::setOrientation(Quaternion orientation) {
    m_orientation = orientation;
}

void StateVector::setAngVelocity(Vector3<dimension::AngularVelocity> ang_velocity) {
    try {
        if (ang_velocity.frame() != ReferenceFrame::BODY)
            throw std::invalid_argument("StateVector: The angular velocity must be specified in the body-fixed frame.");

        m_ang_velocity = ang_velocity;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    m_ang_velocity = ang_velocity;
}

void StateVector::setPositionDerivative(Vector3<dimension::Velocity> velocity) {
    m_position = Vector3<dimension::Distance>(velocity.frame(), {velocity[0], velocity[1], velocity[2]});
}

void StateVector::setVelocityDerivative(Vector3<dimension::Acceleration> acceleration) {
    m_velocity = Vector3<dimension::Velocity>(acceleration.frame(), {acceleration[0], acceleration[1], acceleration[2]});
}

void StateVector::setAngVelocityDerivative(Vector3<dimension::AngularAcceleration> ang_acceleration) {
    m_ang_velocity = Vector3<dimension::AngularVelocity>(ang_acceleration.frame(), {ang_acceleration[0], ang_acceleration[1], ang_acceleration[2]});
}

StateVector StateVector::transformToFrame(ReferenceFrame frame) const {
    if (m_frame == frame)
        return *this;
    if (m_frame == ReferenceFrame::NONE || frame == ReferenceFrame::NONE)
        return *this;

    // Transform the position and velocity vectors
    std::string in_frame(m_frame);
    std::string out_frame(frame);
    double rot[6][6];
    sxform_c(in_frame.c_str(), out_frame.c_str(), m_ephemeris_time, rot);

    double in[6] = {m_position[0], m_position[1], m_position[2], m_velocity[0], m_velocity[1], m_velocity[2]};
    double out[6];
    mxvg_c(rot, in, 6, 6, out);

    Vector3<dimension::Distance> pos_out(frame, {out[0], out[1], out[2]});
    Vector3<dimension::Velocity> vel_out(frame, {dimension::Velocity(out[3]), dimension::Velocity(out[4]), dimension::Velocity(out[5])});
    return StateVector(m_ephemeris_time, pos_out, vel_out, m_orientation, m_ang_velocity);
}

Vector3d StateVector::computeGeodeticCoord() const {
    StateVector itrf93_state = this->transformToFrame(ReferenceFrame::ITRF93);
    Vector3<dimension::Distance> itrf93_position = itrf93_state.position();
    SpiceDouble lon = 0.0, lat = 0.0, alt = 0.0;
    SpiceDouble ecef_position[3] = {itrf93_position[0], itrf93_position[1], itrf93_position[2]};
    recgeo_c(ecef_position, constants::R_EARTH, constants::EARTH_FLATTENING, &lon, &lat, &alt);
    return Vector3d({lat, lon, alt});
}

double StateVector::normInf() const {
    return std::max(std::max(m_position.normInf(), m_velocity.normInf()), std::max(m_orientation.normInf(), m_ang_velocity.normInf()));
}

void StateVector::normalizeQuaternion() {
    m_orientation.normalize();
}


double StateVector::operator[](std::size_t i) const {
    if (i < 3)
        return m_position[i];
    else if (i < 6)
        return m_velocity[i-3];
    else if (i < 10)
        return m_orientation[i-6];
    else
        return m_ang_velocity[i-10];
}

StateVector& StateVector::operator+=(const StateVector& state) {
    m_position += state.position();
    m_velocity += state.velocity();
    m_orientation += state.orientation();
    m_ang_velocity += state.angVelocity();
    return *this;
}

StateVector& StateVector::operator*=(double a) {
    m_position *= a;
    m_velocity *= a;
    m_orientation *= a;
    m_ang_velocity *= a;
    return *this;
}


StateVector operator+(const StateVector& state1, const StateVector& state2) {
    return StateVector(state1.time(), state1.position()+state2.position(), state1.velocity()+state2.velocity(), state1.orientation()+state2.orientation(),
                       state1.angVelocity() +
                                                                                                                                                           state2.angVelocity());
}

StateVector operator*(double a, const StateVector& state) {
    return StateVector(state.time(), a*state.position(), a*state.velocity(), a*state.orientation(), a*
                                                                                                    state.angVelocity());
}

StateVector operator/(const StateVector& state1, const StateVector& state2) {
    return StateVector(state1.time(), state1.position()/state2.position(), state1.velocity()/state2.velocity(), state1.orientation()/state2.orientation(),
                       state1.angVelocity() /
                                                                                                                                                           state2.angVelocity());
}

StateVector abs(const StateVector &state) {
    return StateVector(state.time(), state.position().abs(), state.velocity().abs(), state.orientation().abs(),
                       state.angVelocity().abs());
}

std::ofstream& operator<<(std::ofstream &out, const StateVector& state) {
    std::stringstream ss;
    ss << state.m_ephemeris_time;
    ss << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < 13; ++i) {
        ss << ", ";
        ss << state[i];
    }
    out << ss.rdbuf();
    return out;
}