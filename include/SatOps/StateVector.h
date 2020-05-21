//
// Created by Sylvain  on 4/12/20.
//

#ifndef SATOPS_STATEVECTOR_H
#define SATOPS_STATEVECTOR_H

#include <fstream>
#include <string>
#include <boost/numeric/odeint.hpp>
#include "ReferenceFrame.h"
#include "dimensions.h"
#include "Vector.h"
#include "Quaternion.h"

/** State vector of a spacecraft containing its position, velocity, orientation, and angular velocity.
 *
 * The state vector is composed of 13 elements: 3 for the position, 3 for the velocity, 4 for the orientation (quaternion), and 3 for the angular velocity.
 */
class StateVector {
public:
    StateVector();

    /**
     * @param ephemeris_time    Ephemeris time in seconds.
     * @param position          Vector containing the position in km.
     * @param velocity          Vector containing the velocity in km/s.
     * @param orientation       Quaternion describing the rotation from a base frame to the body-fixed frame.
     * @param ang_velocity      Vector containing the angular velocity in rad/s.
     */
    StateVector(double ephemeris_time, Vector3<dimension::Distance> position, Vector3<dimension::Velocity> velocity, Quaternion orientation, Vector3<dimension::AngularVelocity> ang_velocity);

    /**
     * @param epoch             Initial epoch with format "YYYY-MM-DD hh:mm:ss".
     * @param position          Vector containing the position in km.
     * @param velocity          Vector containing the velocity in km/s.
     * @param orientation       Quaternion describing the rotation from a base frame to the body-fixed frame.
     * @param ang_velocity      Vector containing the angular velocity in rad/s.
     */
    StateVector(const std::string& epoch, Vector3<dimension::Distance> position, Vector3<dimension::Velocity> velocity, Quaternion orientation, Vector3<dimension::AngularVelocity> ang_velocity);


    [[nodiscard]] ReferenceFrame frame() const;
    [[nodiscard]] double time() const;
    [[nodiscard]] Vector3<dimension::Distance> position() const;
    [[nodiscard]] Vector3<dimension::Velocity> velocity() const;
    [[nodiscard]] Quaternion orientation() const;
    [[nodiscard]] Vector3<dimension::AngularVelocity> angVelocity() const;

    void setFrame(ReferenceFrame frame);
    void setTime(double time);
    void setPosition(Vector3<dimension::Distance> position);
    void setVelocity(Vector3<dimension::Velocity> velocity);
    void setOrientation(Quaternion orientation);
    void setAngVelocity(Vector3<dimension::AngularVelocity> ang_velocity);

    void setPositionDerivative(Vector3<dimension::Velocity> velocity);
    void setVelocityDerivative(Vector3<dimension::Acceleration> acceleration);
    void setAngVelocityDerivative(Vector3<dimension::AngularAcceleration> ang_acceleration);

    /** Transforms the position and velocity vector into a given reference frame.
     *
     * @param frame ReferenceFrame to transform into.
     * @return StateVector where the position and velocity are expressed in the given frame.
     */
    [[nodiscard]] StateVector transformToFrame(ReferenceFrame frame) const;

    /** Computes the position in geodetic coordinates.
     *
     * @return Vector containing the geodetic latitude (rad), geodetic longitude (rad), and altitude above the reference spheroid (km).
     */
    [[nodiscard]] Vector3d computeGeodeticCoord() const;

    /** Computes the infinity norm of the state vector. */
    [[nodiscard]] double normInf() const;

    /** Normalizes the quaternion describing the orientation of the spacecraft. */
    void normalizeQuaternion();

    double operator[](std::size_t i) const;
    StateVector& operator+=(const StateVector& state);
    StateVector& operator*=(double a);
    friend StateVector operator+(const StateVector& state1, const StateVector& state2);
    friend StateVector operator*(double a, const StateVector& state);
    friend StateVector operator/(const StateVector& state1, const StateVector& state2);
    friend StateVector abs(const StateVector &state);
    friend std::ofstream& operator<<(std::ofstream &out, const StateVector& state);

private:
    ReferenceFrame m_frame;
    double m_ephemeris_time;
    Vector3<dimension::Distance> m_position;
    Vector3<dimension::Velocity> m_velocity;
    Quaternion m_orientation;
    Vector3<dimension::AngularVelocity> m_ang_velocity;
};


// Template specialization to be compatible with boost odeint
namespace boost::numeric::odeint {
    template<>
    struct vector_space_norm_inf< StateVector >
    {
        typedef double result_type;
        double operator()( const StateVector &state ) const
        {
            return state.normInf();
        }
    };
}

#endif //SATOPS_STATEVECTOR_H
