//
// Created by Sylvain Renevey on 1/26/18.
//



#include "SatOps/EnvironmentModel.h"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <boost/math/special_functions/factorials.hpp>
#include "SatOps/constants.h"
#include "SatOps/StateVector.h"
#include "transformations.h"
extern "C" {
#include "SpiceUsr.h"
};


EnvironmentModel::EnvironmentModel():
        m_central_body(Body::Earth),
        m_gp_degree(0),
        m_c_geopot_coeffs(nullptr),
        m_s_geopot_coeffs(nullptr),
        m_atm_model(AtmModel::None),
        m_third_bodies({}),
        m_earth_gram_atm_model(nullptr),
        m_srp_flag(false),
        m_mag_flag(false),
        m_g_mag_coeffs(nullptr),
        m_h_mag_coeffs(nullptr)
{
}


EnvironmentModel::EnvironmentModel(Body central_body, int gp_degree, const std::string& geopot_model_path, const std::vector<Body>& third_bodies, AtmModel atm_model, const std::string& earthgram_path, const std::string& epoch, bool srp_flag, bool mag_flag, const std::string& igrf_model_path):
        m_central_body(central_body),
        m_gp_degree(gp_degree),
        m_c_geopot_coeffs(nullptr),
        m_s_geopot_coeffs(nullptr),
        m_third_bodies(third_bodies),
        m_atm_model(atm_model),
        m_earth_gram_atm_model(nullptr),
        m_srp_flag(srp_flag),
        m_mag_flag(mag_flag),
        m_g_mag_coeffs(nullptr),
        m_h_mag_coeffs(nullptr)
{
    // Load geopotential coefficients from model file
    if (gp_degree > 1 && !geopot_model_path.empty())
        loadGeopotCoeff(gp_degree, geopot_model_path);

    // Load geomagnetic coefficients
    if (mag_flag)
        loadGeomagCoeff(igrf_model_path);

    // If an atmospheric model is selected, initialize coefficients
    if (m_atm_model == AtmModel::EarthGRAM && m_central_body == Body::Earth) {
        if (!epoch.empty()) {
            m_earth_gram_atm_model = std::make_unique<Atm1>();
            m_earth_gram_atm_model->initdata(earthgram_path, "NameRef.txt", epoch);
        } else {
            std::cerr << "The epoch must be passed in as argument." << std::endl;
        }
    } else if (m_atm_model == AtmModel::Exponential && m_central_body == Body::Earth) {
        initExpModel();
    } else if (m_atm_model != AtmModel::None && m_central_body != Body::Earth) {
        m_atm_model = AtmModel::None;
        std::cerr << "The drag effects are currently only supported for the Earth." << std::endl;
    }
}


Body EnvironmentModel::centralBody() const {
    return m_central_body;
}

std::vector<Body> EnvironmentModel::thirdBodies() const {
    return m_third_bodies;
}

int EnvironmentModel::gpDegree() const {
    return m_gp_degree;
}

int EnvironmentModel::getIndex(int degree, int order) const {
    return degree * (degree + 1) / 2 + order;
}

double EnvironmentModel::cCoeff(int degree, int order) const {
    int i = getIndex(degree, order);
    return m_c_geopot_coeffs[i];
}

double EnvironmentModel::sCoeff(int degree, int order) const {
    int i = getIndex(degree, order);
    return m_s_geopot_coeffs[i];
}


bool EnvironmentModel::isDrag() const {
    return m_atm_model != AtmModel::None;
}


bool EnvironmentModel::isSrp() const {
    return m_srp_flag;
}

bool EnvironmentModel::isMagField() const {
    return m_mag_flag;
}

void EnvironmentModel::loadGeopotCoeff(int degree, std::string model_file_name) {

    // Allocates arrays containing the coefficients
    int num_entries = (degree + 1) * (degree + 2) / 2;
    m_c_geopot_coeffs = std::make_unique<double[]>(num_entries);
    m_s_geopot_coeffs = std::make_unique<double[]>(num_entries);

    for (std::size_t i = 0; i < 3; ++i) {
        m_c_geopot_coeffs[i] = 0.;
        m_s_geopot_coeffs[i] = 0.;
    }

    // Opens file containing the geopotential coefficients and populates arrays
    ifstream data_file(model_file_name);
    if (data_file.is_open()) {

        // Reads required lines
        for (std::size_t i = 3; i < num_entries; ++i) {
            string line;
            getline(data_file, line);

            // Parses the line, converts to double, and stores into arrays
            unsigned long n = line.find('D');
            while (n != string::npos) {
                line.replace(n, 1, "E");
                n = line.find('D');
            }
            std::istringstream iss(line);
            std::vector<std::string> results(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
            m_c_geopot_coeffs[i] = stod(results[2]);
            m_s_geopot_coeffs[i] = stod(results[3]);
        }
        data_file.close();
    } else {
    	std::cerr << "Error opening the file containing the geopotential coefficients." << std::endl;
    }
}

void EnvironmentModel::loadGeomagCoeff(std::string model_file_name) {
    // Allocates arrays containing the coefficients
    int degree = 13;
    int num_entries = (degree + 1) * (degree + 2) / 2;
    m_g_mag_coeffs = std::make_unique<double[]>(num_entries);
    m_h_mag_coeffs = std::make_unique<double[]>(num_entries);

    for (std::size_t i = 0; i < num_entries; ++i) {
        m_g_mag_coeffs[i] = 0.;
        m_h_mag_coeffs[i] = 0.;
    }

    // Opens file containing the geomagnetic coefficients and populates arrays
    ifstream data_file(model_file_name);
    if (data_file.is_open()) {

        // Skip first four lines
        std::string str;
        for (int i = 0; i < 4; ++i) {
            std::getline(data_file, str);
        }

        // Reads required lines
        while (std::getline(data_file, str))
        {
            std::istringstream iss(str);
            std::string token;
            int count = 0;
            std::string coeff = "";
            int degree = 0;
            int order = 0;
            double value = 0.;
            while (iss >> token)
            {
                if (count == 0)
                    coeff = token;
                else if (count == 1)
                    degree = stoi(token);
                else if (count == 2)
                    order = stoi(token);
                else if (count == 27)
                    value = stod(token);

                ++count;
            }

            int i = getIndex(degree, order);
            if (coeff == "g")
                m_g_mag_coeffs[i] = value * 1E-9;
            else if (coeff == "h")
                m_h_mag_coeffs[i] = value * 1E-9;
        }
        data_file.close();
    } else {
        std::cerr << "Error opening the file containing the geomagnetic coefficients." << std::endl;
    }
}

double EnvironmentModel::atmDensity(const StateVector& state, double elapsed_time, double et) const {

    double density = 0.0;
    SpiceDouble inertial_position[3];
    for (std::size_t i = 0; i < 3; ++i)
        inertial_position[i] = state.position()[i];

    if (m_atm_model == AtmModel::EarthGRAM) {

        // Transform the state to ITRF93 (rotating)
        StateVector itrf93_state = state.transformToFrame(ReferenceFrame::ITRF93);
        Vector3<dimension::Distance> itrf93_position = itrf93_state.position();

		// Transform from rectangular to geodetic coordinates
		SpiceDouble lon = 0.0, lat = 0.0, alt = 0.0;
		SpiceDouble ecef_position[3] = {itrf93_position[0], itrf93_position[1], itrf93_position[2]};
		recgeo_c(ecef_position, constants::R_EARTH, constants::EARTH_FLATTENING, &lon, &lat, &alt);
		lat = lat * constants::RAD_TO_DEG;
        lon = lon * constants::RAD_TO_DEG;

        double pm1, dm1, tm1, um1, vm1, wm1, pp1, dp1, tp1, up1, vp1, wp1, ps1, ds1, ts1,
                us1, vs1, ws1, psmall, dsmall, tsmall, usmall, vsmall, wsmall, sos, sosp;

        int iupdate = 1;
        int initonce = 0;

        if (elapsed_time == 0.0)
            initonce = 1;

        m_earth_gram_atm_model->traj(alt, lat, lon, elapsed_time, iupdate, initonce, &dm1, &pm1, &tm1, &um1, &vm1, &wm1,
                                     &dp1, &pp1, &tp1, &up1, &vp1, &wp1, &ds1, &ps1, &ts1, &us1, &vs1, &ws1, &dsmall,
                                     &psmall, &tsmall, &usmall, &vsmall, &wsmall, &sos, &sosp);
        density = dm1;

    } else if (m_atm_model == AtmModel::Exponential) {

		// Compute the altitude above the ellipsoid. There is no atmosphere above 1100 km.
		double h_ell = state.position().norm() - constants::R_EARTH;
        if (h_ell <= 1100.0) {
            for (std::size_t i = 0; i < 28; ++i) {
                if (m_exp_atm_model[i][0] <= h_ell && m_exp_atm_model[i][1] > h_ell) {
                    double h0 = m_exp_atm_model[i][0];
                    double rho0 = m_exp_atm_model[i][2];
                    double H = m_exp_atm_model[i][3];
                    density = rho0 * exp(- (h_ell - h0) / H);
                }
            }
        } else {
            density = 0.0;
        }
    }
    density = density * 1E9; // convert from kg/m^3 to kg/km^3
    return density;
}

Vector3<dimension::Distance> EnvironmentModel::bodyVector(const Body& body, double et, ReferenceFrame frame) const {
    std::string central_body = std::to_string(m_central_body.naifId());
    std::string third_body = std::to_string(body.naifId());

    SpiceDouble third_body_state[6];
    SpiceDouble lt = 0;
    spkezr_c(third_body.c_str(), et, std::string(frame).c_str(), "NONE", central_body.c_str(), third_body_state, &lt);
    return Vector3<dimension::Distance>(frame, third_body_state);
}


std::vector<std::array<double, 2>> EnvironmentModel::geopotHarmonicCoeff(const StateVector& state, double et) const {

    if (state.time() != et)
        throw;

    // Transforms the position to ITRF93 (rotating).
    StateVector state_itrf93 = state.transformToFrame(ReferenceFrame::ITRF93);
    Vector3<dimension::Distance> ecef_position = state_itrf93.position();

    // Computes V and W coefficients
    int num_rows = (m_gp_degree + 2) * (m_gp_degree + 3) / 2;
    double r = ecef_position.norm();

    std::vector<std::array<double, 2>> vw_coeffs(num_rows);
    vw_coeffs[0][0] = constants::R_EARTH_EGM08 / r;
    vw_coeffs[0][1] = 0.;

    vw_coeffs[getIndex(1, 0)][0] = ecef_position[2] * constants::R_EARTH_EGM08 / (r * r) * vw_coeffs[getIndex(0, 0)][0];
    vw_coeffs[getIndex(1, 0)][1] = 0;

    for (int n = 2; n < m_gp_degree + 2; ++n) {
        vw_coeffs[getIndex(n, 0)][0] =
                (2 * n-1) * ecef_position[2] * constants::R_EARTH_EGM08 / (n * r * r) * vw_coeffs[getIndex(n - 1, 0)][0] -
                (n-1) * constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08 / (n * r * r) * vw_coeffs[getIndex(n - 2, 0)][0];
        vw_coeffs[getIndex(n, 0)][1] = 0;
    }

    // Tesseral terms
    for (int n = 1; n < m_gp_degree + 2; ++n) {
        vw_coeffs[getIndex(n, n)][0] = (2 * n - 1) *
                                       (ecef_position[0] * constants::R_EARTH_EGM08 / (r * r) * vw_coeffs[getIndex(n - 1,
                                                                                                                n - 1)][0] -
                                     ecef_position[1] * constants::R_EARTH_EGM08 / (r * r) * vw_coeffs[getIndex(n - 1,
                                                                                                                n - 1)][1]);
        vw_coeffs[getIndex(n, n)][1] = (2 * n - 1) *
                                       (ecef_position[0] * constants::R_EARTH_EGM08 / (r * r) * vw_coeffs[getIndex(n - 1,
                                                                                                                n - 1)][1] +
                                     ecef_position[1] * constants::R_EARTH_EGM08 / (r * r) * vw_coeffs[getIndex(n - 1,
                                                                                                                n - 1)][0]);
    }

    // Remaining terms
    for (int n = 1; n < m_gp_degree + 2; ++n) {
        for (int m = 1; m < n; ++m) {
            if (n == m + 1) {
                vw_coeffs[getIndex(n, m)][0] = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                               vw_coeffs[getIndex(n - 1, m)][0];
                vw_coeffs[getIndex(n, m)][1] = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                               vw_coeffs[getIndex(n - 1, m)][1];
            } else {
                vw_coeffs[getIndex(n, m)][0] = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                               vw_coeffs[getIndex(n - 1, m)][0] -
                                            (n+m-1) * constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08 / ((n-m) * r * r) *
                                            vw_coeffs[getIndex(n - 2, m)][0];
                vw_coeffs[getIndex(n, m)][1] = (2 * n - 1) * ecef_position[2] * constants::R_EARTH_EGM08 / ((n - m) * r * r) *
                                               vw_coeffs[getIndex(n - 1, m)][1] -
                                            (n+m-1) * constants::R_EARTH_EGM08 * constants::R_EARTH_EGM08 / ((n-m) * r * r) *
                                            vw_coeffs[getIndex(n - 2, m)][1];
            }
        }
    }
    return vw_coeffs;
}


Vector3<dimension::Distance> EnvironmentModel::sunSpacecraftVector(const StateVector& state, double et) const {
    Vector3<dimension::Distance> r_sun2sc(state.frame());

    if (m_central_body == Body::Sun) {
        r_sun2sc = state.position();
    }
    else {
        std::string cb = std::to_string(m_central_body.naifId());

        // Retrieve the position vector from the Sun to the central body
        SpiceDouble central_body_state[6];
        SpiceDouble lt = 0.0;
        spkezr_c(cb.c_str(), et, std::string(state.frame()).c_str(), "NONE", "Sun", central_body_state, &lt);
        Vector3<dimension::Distance> r_sun2cb(state.frame(), central_body_state);

        r_sun2sc = r_sun2cb + state.position();
    }
    return r_sun2sc;
}

double EnvironmentModel::isInShadow(const StateVector& state, double et) const {
    Vector3<dimension::Distance> r_sun2sc = sunSpacecraftVector(state, et);

    double a = asin(constants::R_SUN / r_sun2sc.norm());
    double b = asin(m_central_body.radius() / state.position().norm());
    double c = acos( state.position().dot(r_sun2sc) / (state.position().norm() * r_sun2sc.norm()));

    double nu = 0.0;
    if (abs(a-b) < c && c < a + b) {
        double x = (c * c + a * a - b * b) / (2.0 * c);
        double y = sqrt(a * a - x * x);
        double A = a * a * acos(x / a) + b * b * acos((c - x) / b) - c * y;
        nu = 1.0 - A / (constants::PI * a * a);
    } else if (c >= a + b) {
        nu = 0.0;
    } else if (c < b - a) {
        nu = 1.0;
    }
    return nu;
}

void EnvironmentModel::initExpModel() {
    // Second dimension is:
    // - h_ellp (km) lower bound (which corresponds to base altitude h0)
    // - h_ellp (km) upper bound
    // - nominal density (kg/m^3)
    // - scale height H (km)

    m_exp_atm_model[0][0] = 0.0;
    m_exp_atm_model[0][1] = 25.0;
    m_exp_atm_model[0][2] = 1.225;
    m_exp_atm_model[0][3] = 7.249;

    m_exp_atm_model[1][0] = 25.0;
    m_exp_atm_model[1][1] = 30.0;
    m_exp_atm_model[1][2] = 3.899E-2;
    m_exp_atm_model[1][3] = 6.349;

    m_exp_atm_model[2][0] = 30.0;
    m_exp_atm_model[2][1] = 40.0;
    m_exp_atm_model[2][2] = 1.774E-2;
    m_exp_atm_model[2][3] = 6.682;

    m_exp_atm_model[3][0] = 40.0;
    m_exp_atm_model[3][1] = 50.0;
    m_exp_atm_model[3][2] = 3.972E-3;
    m_exp_atm_model[3][3] = 7.554;

    m_exp_atm_model[4][0] = 50.0;
    m_exp_atm_model[4][1] = 60.0;
    m_exp_atm_model[4][2] = 1.057E-3;
    m_exp_atm_model[4][3] = 8.382;

    m_exp_atm_model[5][0] = 60.0;
    m_exp_atm_model[5][1] = 70.0;
    m_exp_atm_model[5][2] = 3.206E-4;
    m_exp_atm_model[5][3] = 7.714;

    m_exp_atm_model[6][0] = 70.0;
    m_exp_atm_model[6][1] = 80.0;
    m_exp_atm_model[6][2] = 8.770E-5;
    m_exp_atm_model[6][3] = 6.549;

    m_exp_atm_model[7][0] = 80.0;
    m_exp_atm_model[7][1] = 90.0;
    m_exp_atm_model[7][2] = 1.905E-5;
    m_exp_atm_model[7][3] = 5.799;

    m_exp_atm_model[8][0] = 90.0;
    m_exp_atm_model[8][1] = 100.0;
    m_exp_atm_model[8][2] = 3.396E-6;
    m_exp_atm_model[8][3] = 5.382;

    m_exp_atm_model[9][0] = 100.0;
    m_exp_atm_model[9][1] = 110.0;
    m_exp_atm_model[9][2] = 5.297E-7;
    m_exp_atm_model[9][3] = 5.877;

    m_exp_atm_model[10][0] = 110.0;
    m_exp_atm_model[10][1] = 120.0;
    m_exp_atm_model[10][2] = 9.661E-8;
    m_exp_atm_model[10][3] = 7.263;

    m_exp_atm_model[11][0] = 120.0;
    m_exp_atm_model[11][1] = 130.0;
    m_exp_atm_model[11][2] = 2.438E-8;
    m_exp_atm_model[11][3] = 9.473;

    m_exp_atm_model[12][0] = 130.0;
    m_exp_atm_model[12][1] = 140.0;
    m_exp_atm_model[12][2] = 8.484E-9;
    m_exp_atm_model[12][3] = 12.636;

    m_exp_atm_model[13][0] = 140.0;
    m_exp_atm_model[13][1] = 150.0;
    m_exp_atm_model[13][2] = 3.845E-9;
    m_exp_atm_model[13][3] = 16.149;

    m_exp_atm_model[14][0] = 150.0;
    m_exp_atm_model[14][1] = 180.0;
    m_exp_atm_model[14][2] = 2.070E-9;
    m_exp_atm_model[14][3] = 22.523;

    m_exp_atm_model[15][0] = 180.0;
    m_exp_atm_model[15][1] = 200.0;
    m_exp_atm_model[15][2] = 5.464E-10;
    m_exp_atm_model[15][3] = 29.740;

    m_exp_atm_model[16][0] = 200.0;
    m_exp_atm_model[16][1] = 250.0;
    m_exp_atm_model[16][2] = 2.789E-10;
    m_exp_atm_model[16][3] = 37.105;

    m_exp_atm_model[17][0] = 250.0;
    m_exp_atm_model[17][1] = 300.0;
    m_exp_atm_model[17][2] = 7.248E-11;
    m_exp_atm_model[17][3] = 45.546;

    m_exp_atm_model[18][0] = 300.0;
    m_exp_atm_model[18][1] = 350.0;
    m_exp_atm_model[18][2] = 2.418E-11;
    m_exp_atm_model[18][3] = 53.628;

    m_exp_atm_model[19][0] = 350.0;
    m_exp_atm_model[19][1] = 400.0;
    m_exp_atm_model[19][2] = 9.518E-12;
    m_exp_atm_model[19][3] = 53.298;

    m_exp_atm_model[20][0] = 400.0;
    m_exp_atm_model[20][1] = 450.0;
    m_exp_atm_model[20][2] = 3.725E-12;
    m_exp_atm_model[20][3] = 58.515;

    m_exp_atm_model[21][0] = 450.0;
    m_exp_atm_model[21][1] = 500.0;
    m_exp_atm_model[21][2] = 1.585E-12;
    m_exp_atm_model[21][3] = 60.828;

    m_exp_atm_model[22][0] = 500.0;
    m_exp_atm_model[22][1] = 600.0;
    m_exp_atm_model[22][2] = 6.967E-13;
    m_exp_atm_model[22][3] = 63.822;

    m_exp_atm_model[23][0] = 600.0;
    m_exp_atm_model[23][1] = 700.0;
    m_exp_atm_model[23][2] = 1.454E-13;
    m_exp_atm_model[23][3] = 71.835;

    m_exp_atm_model[24][0] = 700.0;
    m_exp_atm_model[24][1] = 800.0;
    m_exp_atm_model[24][2] = 3.614E-14;
    m_exp_atm_model[24][3] = 88.667;

    m_exp_atm_model[25][0] = 800.0;
    m_exp_atm_model[25][1] = 900.0;
    m_exp_atm_model[25][2] = 1.170E-14;
    m_exp_atm_model[25][3] = 124.64;

    m_exp_atm_model[26][0] = 900.0;
    m_exp_atm_model[26][1] = 1000.0;
    m_exp_atm_model[26][2] = 5.245E-15;
    m_exp_atm_model[26][3] = 181.05;

    m_exp_atm_model[27][0] = 1000.0;
    m_exp_atm_model[27][1] = 1100.0;
    m_exp_atm_model[27][2] = 3.019E-15;
    m_exp_atm_model[27][3] = 268.0;
}

Vector3d EnvironmentModel::magneticField(const StateVector& state, double et) const {

    // Creates lambdas to retrieve C and S coefficients from the environment model
    auto C = [&](int degree, int order) { return m_g_mag_coeffs[getIndex(degree, order)]; };
    auto S = [&](int degree, int order) { return m_h_mag_coeffs[getIndex(degree, order)]; };

    // Computes and retrieves the V and W coefficients
    auto vw_coeffs = geopotHarmonicCoeff(state, et);
    auto V = [&](int degree, int order) { return vw_coeffs[getIndex(degree, order)][0]; };
    auto W = [&](int degree, int order) { return vw_coeffs[getIndex(degree, order)][1]; };


    Vector3d B(ReferenceFrame::ITRF93);

    // n = degree, m = order
    for (int n = 0; n < 13 + 1; ++n) {
        for (int m = 0; m < n + 1; ++m) {
            double scale = 0;
            if (m == 0) {
                scale = sqrt((2 * n + 1));
                B[0] += scale * (-C(n, 0) * V(n+1, 1));
                B[1] += scale * (-C(n, 0) * W(n+1, 1));
            } else {
                scale = sqrt(2 * (2 * n + 1) * boost::math::factorial<double>(n-m) / boost::math::factorial<double>(n+m));
                B[0] += scale * 0.5 * ((-C(n, m) * V(n+1, m+1) - S(n, m) * W(n+1, m+1)) +
                                boost::math::factorial<double>(n-m+2) / boost::math::factorial<double>(n-m) *
                                (C(n, m) * V(n+1, m-1) + S(n, m) * W(n+1, m-1)));
                B[1] += scale * 0.5 * ((-C(n, m) * W(n+1, m+1) + S(n, m) * V(n+1, m+1)) +
                                boost::math::factorial<double>(n-m+2) / boost::math::factorial<double>(n-m) *
                                (-C(n, m) * W(n+1, m-1) + S(n, m) * V(n+1, m-1)));
            }
            B[2] += scale * ((n-m+1) * (-C(n, m) * V(n+1, m) - S(n, m) * W(n+1, m)));
        }
    }
    return transformations::rotateToFrame(B, state.frame(), et);

    /*
    Vector3d<double> m(pow(constants::IGRF_2020_RADIUS, 3) * Vector3d<double>(ReferenceFrame::ITRF93, {constants::IGRF_2020_G11, constants::IGRF_2020_H11, constants::IGRF_2020_G01}));

    // Transform the position to ITRF93 (rotating).
    StateVector state_itrf93 = state.transformToFrame(ReferenceFrame::ITRF93);
    Vector3d<Dimension::Distance> r_itrf93 = state_itrf93.position();

    // Compute the magnetic field
    Vector3d<double> B = 3.0 * (m.dot(r_itrf93) * r_itrf93 - r_itrf93.norm()*r_itrf93.norm()*m) / pow(r_itrf93.norm(), 5);

    // Transform back to initial frame and return
    return transformations::rotateToFrame(B, state.frame(), et);
     */
}