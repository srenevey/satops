//
// Created by Sylvain  on 4/7/20.
//

#include "SatOps/Magnetometer.h"
#include "transformations.h"
#include "SatOps/EnvironmentModel.h"
#include "SatOps/StateVector.h"
#include <cstddef>


Magnetometer::Magnetometer(const EnvironmentModel& env_model, Vector3d bias, Vector3d stddev, Matrix3d rotation_matrix, double max_value, double min_value, double resolution): m_env_model(env_model), m_rot_bff2sff(rotation_matrix), m_max_value(max_value), m_min_value(min_value), m_resolution(resolution) {
    for (std::size_t i = 0; i < 3; ++i) {
        m_distribution[i] = std::normal_distribution<double>(bias[i], stddev[i]);
    }
}

Magnetometer::~Magnetometer() {}

Vector3d Magnetometer::measure(const StateVector& state) {
    Vector3d B = m_env_model.magneticField(state, state.time());

    Vector3d noise(state.frame(), {m_distribution[0](m_generator), m_distribution[1](m_generator), m_distribution[2](m_generator)});
    Vector3d B_meas = B + noise;
    B_meas.clamp(m_min_value, m_max_value);

    // Quantizes values to resolution level
    B_meas.roundDown(m_resolution);

    std::string filename = "magnetic_data.dat";
    {
        std::ofstream f(filename, std::ios::app);
        f << state.time() << ", " << B[0] << ", " << B[1] << ", " << B[2] << ", " << B_meas[0] << ", " << B_meas[1] << ", " << B_meas[2] << '\n';
    }


    Vector3d B_meas_bff = transformations::rotateToFrame(B_meas, ReferenceFrame::BODY, state.time(),
                                                         state.orientation());
    Vector3d B_meas_sff = m_rot_bff2sff * B_meas_bff;
    return B_meas_sff;
}

Vector3d Magnetometer::bias() const {
    return Vector3d({m_distribution[0].mean(), m_distribution[1].mean(), m_distribution[2].mean()});
}

Vector3d Magnetometer::stddev() const {
    return Vector3d({m_distribution[0].stddev(), m_distribution[1].stddev(), m_distribution[2].stddev()});
}