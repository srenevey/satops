//
// Created by Sylvain  on 5/26/20.
//

#include "SatOps/Gyroscope.h"
#include "SatOps/StateVector.h"
#include "SatOps/ReferenceFrame.h"
#include <iostream>

Gyroscope::Gyroscope(Matrix3d rotation_matrix, Matrix3d covariance_matrix, std::function<Vector3d(double)> drift_rate)
    :
    Sensor(rotation_matrix, covariance_matrix),
    m_covariance_matrix(covariance_matrix),
    m_drift_rate(drift_rate)
{}

Gyroscope::~Gyroscope() {}

Vector3d Gyroscope::measure(const StateVector& state) {
    Vector3d bias_sff({0., 0., 0.});
    if (m_drift_rate) {
        bias_sff = m_drift_rate(state.time());
    }

    Vector3d noise_sff = Vector3d(ReferenceFrame::NONE, {m_distribution[0](m_generator), m_distribution[1](m_generator),
                                   m_distribution[2](m_generator)});

    Vector3d ang_velocity_bff = state.angVelocity();
    Vector3d ang_velocity_meas_sff = m_rot_bff2sff * ang_velocity_bff + bias_sff + noise_sff;

    Vector3d ang_velocity_sff = m_rot_bff2sff * ang_velocity_bff;
    std::string filename = "gyroscope_data.dat";
    {
        std::ofstream f(filename, std::ios::app);
        f << state.time() << ", " << ang_velocity_sff[0] << ", " << ang_velocity_sff[1] << ", " << ang_velocity_sff[2] << ", " << ang_velocity_meas_sff[0] << ", " << ang_velocity_meas_sff[1]
          << ", " << ang_velocity_meas_sff[2] << '\n';
    }

    return ang_velocity_meas_sff;
}