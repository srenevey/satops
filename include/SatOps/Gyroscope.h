//
// Created by Sylvain  on 5/26/20.
//

#ifndef SATOPS_GYROSCOPE_H
#define SATOPS_GYROSCOPE_H

#include "Sensor.h"
#include "Vector.h"
#include "Matrix.h"
#include <array>
#include <random>
#include <functional>

class StateVector;

/** Model of a gyroscope (work in progress) */
class Gyroscope: public Sensor {
public:
    /**
     * @param[in] rotation_matrix       Rotation matrix to go from body-fixed to sensor-fixed frame
     * @param[in] covariance_matrix     Covariance matrix in the sensor-fixed frame used to generate white noise
     * @param[in] drift_rate            Function returning the biases as a function of the absolute ephemeris time in the sensor-fixed frame
     */
    Gyroscope(Matrix3d rotation_matrix, Matrix3d covariance_matrix, std::function<Vector3d(double)> drift_rate = {});
    ~Gyroscope();

    /** Measures the angular rates.
     *
     * @param[in] state     State vector of the spacecraft
     * @return              Angular rates in the sensor-fixed frame (in rad/s).
     */
    [[nodiscard]] Vector3d measure(const StateVector& state);

private:
    Matrix3d m_covariance_matrix;
    std::function<Vector3d(double)> m_drift_rate;
};


#endif //SATOPS_GYROSCOPE_H
