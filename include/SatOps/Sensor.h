//
// Created by Sylvain  on 4/7/20.
//

#ifndef SATOPS_SENSOR_H
#define SATOPS_SENSOR_H

#include "Vector.h"
#include "Matrix.h"
#include <array>
#include <random>

class StateVector;

/** Base class representing a sensor */
class Sensor {
public:
    Sensor(Matrix3d rot_bff2sff, Matrix3d covariance_matrix): m_rot_bff2sff(rot_bff2sff) {
        for (std::size_t i = 0; i < 3; ++i) {
            m_distribution[i] = std::normal_distribution<double>(0.0, covariance_matrix[i][i]);
        }
    }
    virtual ~Sensor() {}

    /**
     * @param[in] state     State vector of the spacecraft.
     * @return              Measurement in the sensor-fixed frame.
     */
    virtual Vector3d measure(const StateVector& state) = 0;

    /** Returns the rotation matrix to go from body-fixed frame to sensor-fixed frame */
    [[nodiscard]] Matrix3d getRotationMatrixBFFToSFF() const { return m_rot_bff2sff; }

protected:
    Matrix3d m_rot_bff2sff; // rotation matrix from body-fixed frame to sensor-fixed frame
    std::default_random_engine m_generator;
    std::array<std::normal_distribution<double>, 3> m_distribution; // normal probability distributions used to generate white noise
};


#endif //SATOPS_SENSOR_H
