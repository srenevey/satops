//
// Created by Sylvain  on 4/7/20.
//

#include "SatOps/Magnetometer.h"
#include "transformations.h"
#include "SatOps/EnvironmentModel.h"
#include "SatOps/ReferenceFrame.h"
#include "SatOps/StateVector.h"
#include <cstddef>


Magnetometer::Magnetometer(const EnvironmentModel& env_model, Matrix3d covariance_matrix, Matrix3d rotation_matrix, double max_value, double min_value, double resolution, std::function<Vector3d(double)> drift_rate)
    :
    Sensor(rotation_matrix, covariance_matrix),
    m_env_model(env_model),
    m_covariance_matrix(covariance_matrix),
    m_max_value(max_value),
    m_min_value(min_value),
    m_resolution(resolution),
    m_drift_rate(drift_rate)
{}

Magnetometer::~Magnetometer() {}

Vector3d Magnetometer::measure(const StateVector& state) {
    Vector3d B_meas_sff({0., 0., 0.});
    if (m_env_model.isMagField()) {
        Vector3d B = m_env_model.computeMagneticField(state);
        Vector3d B_bff = transformations::rotateToFrame(B, ReferenceFrame::BODY, state.time(), state.orientation());
        Vector3d B_sff = m_rot_bff2sff * B_bff;

        Vector3d noise_sff({m_distribution[0](m_generator), m_distribution[1](m_generator),
                                       m_distribution[2](m_generator)});

        Vector3d bias_sff({0., 0., 0.});
        if (m_drift_rate) {
            bias_sff = m_drift_rate(state.time());
        }

        B_meas_sff = B_sff + bias_sff + noise_sff;
        B_meas_sff.clamp(m_min_value, m_max_value);

        // Quantizes values to resolution level
        B_meas_sff.roundDown(m_resolution);

        std::string filename = "magnetic_data.dat";
        {
            std::ofstream f(filename, std::ios::app);
            f << state.time() << ", " << B_sff[0] << ", " << B_sff[1] << ", " << B_sff[2] << ", " << B_meas_sff[0] << ", " << B_meas_sff[1]
              << ", " << B_meas_sff[2] << '\n';
        }
    }
    return B_meas_sff;
}