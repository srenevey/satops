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
#include <functional>

class EnvironmentModel;
class StateVector;

/** Model of a 3-axis magnetometer (work in progress). */
class Magnetometer: public Sensor {
public:
    /**
     * @param[in] env_model             Reference to the environment model
     * @param[in] covariance_matrix     Covariance matrix in the sensor-fixed frame used to generate white noise
     * @param[in] rotation_matrix       Rotation matrix to go from body-fixed to sensor-fixed frame
     * @param[in] max_value             Maximum value that the magnetometer can read in Tesla
     * @param[in] min_value             Minimum value that the magnetometer can read in Tesla
     * @param[in] resolution            Resolution of the magnetometer in Tesla
     * @param[in] drift_rate            Function returning the biases as a function of the absolute ephemeris time in the sensor-fixed frame
     */
    Magnetometer(const EnvironmentModel& env_model, Matrix3d covariance_matrix, Matrix3d rotation_matrix, double max_value, double min_value, double resolution, std::function<Vector3d(double)> drift_rate = {});
    ~Magnetometer();

    /** Measures the local magnetic field.
     *
     * @param[in] state     State vector of the spacecraft
     * @return              Magnetic field in the sensor-fixed frame (in Tesla).
     */
    [[nodiscard]] Vector3d measure(const StateVector& state);

private:
    const EnvironmentModel& m_env_model;
    Matrix3d m_covariance_matrix;
    double m_max_value;
    double m_min_value;
    double m_resolution;
    std::function<Vector3d(double)> m_drift_rate;
};


#endif //SATOPS_MAGNETOMETER_H
