//
// Created by Sylvain  on 6/1/20.
//

#include "SatOps/GravitationalField.h"
#include "SatOps/constants.h"
#include "SatOps/ReferenceFrame.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <boost/math/special_functions/factorials.hpp>

GravitationalField::GravitationalField(Model model, unsigned int max_degree, const std::filesystem::path& file)
    :
    PotentialField(max_degree),
    m_model(model)
{
    // Allocates and initializes vectors containing the coefficients
    unsigned int num_entries = (max_degree + 1) * (max_degree + 2) / 2;
    m_c_coeffs = std::vector(num_entries, 0.0);
    m_s_coeffs = std::vector(num_entries, 0.0);

    switch (m_model) {
        case Model::EGM2008:
            loadEGM08Coeffs(file);
            m_radius = constants::R_EARTH_EGM08;
            m_scale_factor = -constants::MU_EARTH_EGM08;
            break;
        case Model::GEM10:
            if (max_degree > 5)
                throw std::invalid_argument("GravitationalField: The degree of expansion of the GEM10 model is too high (max 5).");

            loadEGM08Coeffs(file);
            m_radius = constants::R_EARTH_GEM10;
            m_scale_factor = -constants::MU_EARTH_GEM10;
            break;
    }
}

void GravitationalField::loadEGM08Coeffs(const std::filesystem::path& path) {
    std::ifstream data_file(path);
    if (data_file.is_open()) {
        unsigned int num_entries = (m_max_degree + 1) * (m_max_degree + 2) / 2;

        // Reads required lines
        for (std::size_t i = 3; i < num_entries; ++i) {
            std::string line;
            getline(data_file, line);

            // Parses the line, converts to double, and stores into arrays
            unsigned long k = line.find('D');
            while (k != std::string::npos) {
                line.replace(k, 1, "E");
                k = line.find('D');
            }
            std::istringstream iss(line);
            std::vector<std::string> results(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
            int n = stoi(results[0]);
            int m = stoi(results[1]);
            double kron_delta = (0 == m) ? 1. : 0.;
            long double unnormalization_factor = std::sqrt(boost::math::factorial<long double>(n - m) * (2. * n + 1) * (2. - kron_delta) / boost::math::factorial<long double>(n + m));
            m_c_coeffs[i] = unnormalization_factor * stod(results[2]);
            m_s_coeffs[i] = unnormalization_factor * stod(results[3]);
        }
        data_file.close();
    } else {
        throw std::invalid_argument("GravitationalField: The file containing the spherical harmonics coefficients could not be opened.");
    }
}

Vector3d GravitationalField::getAcceleration(const Vector3d& itrf93_position, Matrix3d* jacobian) {
    double x = itrf93_position[0];
    double y = itrf93_position[1];
    double z = itrf93_position[2];

    std::array<double, 12> sums{};
    sums.fill(0.0);
    sums[0] = 1.0;
    sums[4] = 2.0;

    Vector3d gravity_field = computeGradient(x, y, z, sums, jacobian);
    gravity_field.setFrame(ReferenceFrame::ITRF93);
    return gravity_field;
}