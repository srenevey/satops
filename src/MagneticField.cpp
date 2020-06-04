//
// Created by Sylvain  on 5/30/20.
//

#include "SatOps/MagneticField.h"
#include "SatOps/constants.h"
#include "SatOps/ReferenceFrame.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <boost/math/special_functions/factorials.hpp>


MagneticField::MagneticField(Model model, unsigned int max_degree, const std::filesystem::path& file, double year)
    :
    PotentialField(max_degree),
    m_model(model)
{
    // Allocates and initializes vectors containing the coefficients
    unsigned int num_entries = (max_degree + 1) * (max_degree + 2) / 2;
    m_c_coeffs = std::vector(num_entries, 0.0);
    m_s_coeffs = std::vector(num_entries, 0.0);

    switch (m_model) {
        case Model::IGRF13:
            if (m_max_degree > 13)
                throw std::invalid_argument("MagneticField: The degree of expansion is out of range (max 13).");
            if (year < 1900. || year > 2025.)
                throw std::invalid_argument("MagneticField: The year is out of range. The model is valide from 1900 to 2025.");

            loadIGRF13Coeffs(file, year);
            m_radius = constants::IGRF13_RADIUS;
            m_scale_factor = 1E-9 * constants::IGRF13_RADIUS * constants::IGRF13_RADIUS;
            break;
    }
}

// The coefficients of the model are computed for the initial epoch of the simulation
void MagneticField::loadIGRF13Coeffs(const std::filesystem::path& path, double year) {

    // Gets reference epoch of the model and difference with desired time
    int T0 = static_cast<int>(year) - static_cast<int>(std::floor(year)) % 5;
    double dt = year - T0;
    int year_col = (T0 - 1900) / 5 + 3;

    // Opens file containing the geomagnetic coefficients and populates arrays
    std::ifstream data_file(path);
    if (data_file.is_open()) {

        // Skips first four lines (header)
        std::string str;
        for (int i = 0; i < 4; ++i) {
            std::getline(data_file, str);
        }

        // Reads required lines
        while (std::getline(data_file, str)) {
            std::istringstream iss(str);
            std::string value;
            std::string coeff;
            int col = 0, degree = 0, order = 0;
            double main_field_coeff = 0., secular_variation_coeff = 0.;

            while (iss >> value) {
                if (col == 0)
                    coeff = value;
                else if (col == 1)
                    degree = stoi(value);
                else if (col == 2)
                    order = stoi(value);
                else if (col == year_col)
                    main_field_coeff = stod(value);
                else if (col == 28) // if last column, secular variation is directly given
                    secular_variation_coeff = stod(value);
                else if (col == year_col + 1)
                    secular_variation_coeff = (stod(value) - main_field_coeff) / 5.0;

                ++col;
            }

            if (degree > m_max_degree)
                break;

            int i = getIndex(degree, order);
            double kron_delta = (0 == order) ? 1. : 0.;
            double unnormalization_factor = std::sqrt(2.0 * boost::math::factorial<double>(degree-order) / boost::math::factorial<double>(degree+order) - kron_delta);
            if (coeff == "g") {
                m_c_coeffs[i] = unnormalization_factor * (main_field_coeff + dt * secular_variation_coeff);
            } else if (coeff == "h") {
                m_s_coeffs[i] = unnormalization_factor * (main_field_coeff + dt * secular_variation_coeff);
            }
        }
        data_file.close();
    } else {
        throw std::invalid_argument("MagneticField: The file containing the magnetic coefficients could not be opened.");
    }
}

Vector3d MagneticField::getMagneticField(const Vector3d& itrf93_position, Matrix3d* jacobian) {
    double x = itrf93_position[0];
    double y = itrf93_position[1];
    double z = itrf93_position[2];
    double r = itrf93_position.norm();

    std::array<double, 12> sums{};
    sums.fill(0.0);
    sums[0] = m_radius / (r*r) * (3. * (m_c_coeffs[getIndex(1,1)] * x + m_s_coeffs[getIndex(1,1)] * y) + 2. * m_c_coeffs[getIndex(1,0)] * z);
    sums[1] = m_radius / r * m_c_coeffs[getIndex(1,0)];
    sums[2] = m_radius / r * m_c_coeffs[getIndex(1,1)];
    sums[3] = m_radius / r * m_s_coeffs[getIndex(1,1)];
    sums[4] = 2.0;

    Vector3d mag_field = computeGradient(x, y, z, sums, jacobian);
    mag_field.setFrame(ReferenceFrame::ITRF93);
    return mag_field;
}