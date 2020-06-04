//
// Created by Sylvain  on 6/2/20.
//

#include "SatOps/Atmosphere.h"
#include "SatOps/constants.h"
#include "SatOps/StateVector.h"
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <vector>


Atmosphere::Atmosphere(Model model, const std::filesystem::path& file, std::string epoch)
    :
    m_model(model),
    m_earth_gram_model(nullptr),
    m_exp_model_base_alt(nullptr),
    m_exp_model_density(nullptr),
    m_exp_model_scale_height(nullptr)
{
    switch (m_model) {
        case Model::earthGRAM:
            if (epoch.empty())
                throw std::invalid_argument("Atmosphere: An initial epoch must be provided.");
            m_earth_gram_model = std::make_unique<Atm1>();
            m_earth_gram_model->initdata(file.string(), "NameRef.txt", epoch);
            break;
        case Model::exp:
            loadExpModel(file);
            break;
        default:
            break;
    }
}

void Atmosphere::loadExpModel(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (file.is_open()) {
        m_exp_model_base_alt = std::make_unique<double[]>(28);
        m_exp_model_density = std::make_unique<double[]>(28);
        m_exp_model_scale_height = std::make_unique<double[]>(28);

        // Reads required lines
        for (std::size_t i = 0; i < 28; ++i) {
            std::string line;
            getline(file, line);

            // Parses the line, converts to double, and stores into arrays
            std::istringstream iss(line);
            std::vector<std::string> results(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
            m_exp_model_base_alt[i] = stod(results[0]);
            m_exp_model_density[i] = stod(results[2]);
            m_exp_model_scale_height[i] = stod(results[3]);
        }
        file.close();
    } else {
        throw std::invalid_argument("Atmosphere: The file containing the coefficients of the model could not be opened.");
    }
}

double Atmosphere::getDensity(const StateVector& state, double elapsed_time) const {
    double density = 0.0;

    switch (m_model) {
        case Model::earthGRAM: {
            Vector3d geodetic_coord = state.computeGeodeticCoord();
            double lat_deg = geodetic_coord[0] * constants::RAD_TO_DEG;
            double lon_deg = geodetic_coord[1] * constants::RAD_TO_DEG;

            double pm1, tm1, um1, vm1, wm1, pp1, dp1, tp1, up1, vp1, wp1, ps1, ds1, ts1,
                    us1, vs1, ws1, psmall, dsmall, tsmall, usmall, vsmall, wsmall, sos, sosp;

            int iupdate = 1;
            int initonce = (elapsed_time == 0.0) ? 1 : 0;
            m_earth_gram_model->traj(geodetic_coord[2], lat_deg, lon_deg, elapsed_time, iupdate, initonce, &density, &pm1,
                                     &tm1, &um1, &vm1, &wm1,
                                     &dp1, &pp1, &tp1, &up1, &vp1, &wp1, &ds1, &ps1, &ts1, &us1, &vs1, &ws1, &dsmall,
                                     &psmall, &tsmall, &usmall, &vsmall, &wsmall, &sos, &sosp);
            break;
        }
        case Model::exp: {
            // Compute the altitude above the ellipsoid. There is no atmosphere above 1100 km.
            double h_ell = state.position().norm() - constants::R_EARTH;
            if (h_ell >= 0. && h_ell <= 1100.) {
                int index = searchAltitude(h_ell);
                double h0 = m_exp_model_base_alt[index];
                double rho0 = m_exp_model_density[index];
                double H = m_exp_model_scale_height[index];
                density = rho0 * std::exp(-(h_ell - h0) / H);
            }
            break;
        }
        default:
            break;
    }
    return density * 1E9; // converts from kg/m^3 to kg/km^3
}

int Atmosphere::searchAltitude(double alt) const {
    int left = 0;
    int right = 28;

    while (left < right) {
        int m = std::floor((left + right) / 2);
        if (m_exp_model_base_alt[m] <= alt) {
            left = m + 1;
        } else {
            right = m;
        }
    }
    return left-1;
}