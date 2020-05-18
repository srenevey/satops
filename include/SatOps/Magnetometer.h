//
// Created by Sylvain  on 4/7/20.
//

#ifndef SATOPS_MAGNETOMETER_H
#define SATOPS_MAGNETOMETER_H

#include "Sensor.h"
#include "Vector.h"
#include "Matrix.h"
#include <array>
#include <random>

class EnvironmentModel;
class StateVector;

/** Model of a 3-axis magnetometer (work in progress). */
class Magnetometer: public Sensor {
public:
    /**
     * @param bias                  Biases on the three axes
     * @param stddev                Standard deviation on the three axes
     * @param rotation_matrix       Rotation matrix to go from body-fixed to sensor-fixed frame
     * @param max_value             Maximum value that the magnetometer can read
     * @param min_value             Minimum value that the magnetometer can read
     */
    Magnetometer(const EnvironmentModel& env_model, Vector3d bias, Vector3d stddev, Matrix3d rotation_matrix, double max_value, double min_value, double resolution);
    ~Magnetometer();

    /** Measures the local magnetic field.
     *
     * @param state         State vector of the spacecraft
     * @return              Magnetic field in the sensor-fixed frame (in Tesla).
     */
    [[nodiscard]] Vector3d measure(const StateVector& state);
    [[nodiscard]] Vector3d bias() const;
    [[nodiscard]] Vector3d stddev() const;

private:
    const EnvironmentModel& m_env_model;
    Matrix3d m_rot_bff2sff; // rotation matrix from body-fixed frame to sensor-fixed frame
    double m_max_value;
    double m_min_value;
    double m_resolution;
    std::default_random_engine m_generator;
    std::array<std::normal_distribution<double>, 3> m_distribution;
};


#endif //SATOPS_MAGNETOMETER_H
