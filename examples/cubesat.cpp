#include <array>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include "SatOps/SatOps.h"

using namespace unit;

int main() {

    // Define paths to the meta-kernel, geopotential model, geomagnetic model, and EarthGRAM lib
    std::filesystem::path meta_kernel("/Users/sylvain/Developer/satops/assets/kernels.tm");
    std::filesystem::path geopot_model_path("/Users/sylvain/Developer/satops/assets/EGM2008_to2190_TideFree");
    std::filesystem::path geomag_model_path("/Users/sylvain/Developer/satops/assets/igrf13coeffs.txt");
    std::filesystem::path earthgram_path("/Users/sylvain/Developer/satops/extern/earthGRAM2016/");

    // Create a Sim object to load the SPICE kernels
    Sim sim(meta_kernel);

    // Define the initial epoch, final time, and step size.
    std::string epoch("2020-04-05 16:00:00");
    dimension::Time t0 = 0_s;
    dimension::Time t1 = 24_h;
    dimension::Time dt = 5_s;

    // Create the environment model
    GravitationalField g_field(GravitationalField::EGM2008, 10, geopot_model_path);
    MagneticField mag_field(MagneticField::IGRF13, 13, geomag_model_path, 2020.5);
    Atmosphere atm(Atmosphere::earthGRAM, earthgram_path, epoch);
    EnvironmentModel env_model(Body::Earth, &g_field, &mag_field, &atm, {Body::Sun, Body::Moon}, true);


    // Create the initial state of the spacecraft
    Vector3<dimension::Distance> position(ReferenceFrame::J2000, {-7009.22977_km, 0.0_km, 0.0_km});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::J2000, {0.0_kms, -7.0033_kms, 3.2657_kms});
    Quaternion orientation(ReferenceFrame::J2000, 0., 0., 0., 1.);
    Vector3<dimension::AngularVelocity> angular_velocity(ReferenceFrame::BODY, {0.2_degs, 0.0_degs, 0.1_degs});
    StateVector initial_state(epoch, position, velocity, orientation, angular_velocity);

    // Model of a 1U CubeSat
    // Spacecraft are represented as a collection of N flat plates.
    dimension::Mass mass = 1.33_kg;
    dimension::Distance side = 0.1_m;
    double I_ii = mass * side * side / 6.;
    Matrix3d inertia_matrix({{{I_ii, 0., 0.}, {0., I_ii, 0.}, {0., 0., I_ii}}});
    std::vector<Vector3d> face_normals;
    face_normals.push_back(Vector3d(ReferenceFrame::BODY, {1., 0., 0.}));
    face_normals.push_back(Vector3d(ReferenceFrame::BODY, {0., 1., 0.}));
    face_normals.push_back(Vector3d(ReferenceFrame::BODY, {0., 0., 1.}));
    face_normals.push_back(Vector3d(ReferenceFrame::BODY, {-1., 0., 0.}));
    face_normals.push_back(Vector3d(ReferenceFrame::BODY, {0., -1., 0.}));
    face_normals.push_back(Vector3d(ReferenceFrame::BODY, {0., 0., -1.}));

    dimension::Area face_area = side * side;
    std::vector<dimension::Area> face_areas{face_area, face_area, face_area, face_area, face_area, face_area};

    std::vector<Vector3<dimension::Distance>> cop_positions;
    cop_positions.push_back(Vector3<dimension::Distance>(ReferenceFrame::BODY, {side / 2., 0._m, 0._m}));
    cop_positions.push_back(Vector3<dimension::Distance>(ReferenceFrame::BODY, {0._m, side / 2., 0._m}));
    cop_positions.push_back(Vector3<dimension::Distance>(ReferenceFrame::BODY, {0._m, 0._m, side / 2.}));
    cop_positions.push_back(Vector3<dimension::Distance>(ReferenceFrame::BODY, {-side / 2., 0._m, 0._m}));
    cop_positions.push_back(Vector3<dimension::Distance>(ReferenceFrame::BODY, {0._m, -side / 2., 0._m}));
    cop_positions.push_back(Vector3<dimension::Distance>(ReferenceFrame::BODY, {0._m, 0._m, -side / 2.}));

    std::vector<double> specular_reflection_coeff{0.042, 0.042, 0.184, 0.042, 0.042, 0.184}; // Solar panel: 0.042. Gold foil: 0.184.
    std::vector<double> diffuse_reflection_coeff{0.168, 0.168, 0.736, 0.168, 0.168, 0.736}; // Solar panel: 0.168. Gold foil: 0.736.

    // Residual dipole
    dimension::Area area(3.1415 * 7._cm * 7._cm);
    double current = 0.1; // A
    double n_turns = 5;
    Vector3d residual_dipole = n_turns * current * area * Vector3d(ReferenceFrame::BODY, {0.8660, 0.5, 0.});

    // Sensors (experimental)
    std::vector<std::shared_ptr<Sensor>> sensors;
    // Magnetometer
    Matrix3d cov_mat_mag = Matrix3d::diag({1.0E-7, 1.0E-7, 1.0E-7});
    Matrix3d rot_bff2sff({{{0., 0.2588, 0.9659}, {-1., 0., 0.}, {0., -0.9659, 0.2588}}}); // rotation matrix from body-fixed frame to sensor-fixed frame
    double max_value = 1.0E-4; // Tesla
    double min_value = -1.0E-4; // Tesla
    double resolution = 0.2E-9; // Tesla
    auto biases = [](double){ return Vector3d({2.0E-9, 2.0E-9, 2.0E-9}); }; // constant biases in the sensor-fixed frame in Tesla
    std::shared_ptr<Magnetometer> mag1 = std::make_shared<Magnetometer>(env_model, cov_mat_mag, rot_bff2sff, max_value, min_value, resolution, biases);
    sensors.push_back(std::move(mag1));

    // Gyroscope
    Matrix3d cov_mat_gyro = Matrix3d::diag({1.0E-4, 1.0E-4, 1.0E-4});
    rot_bff2sff = Matrix3d({{{0., 0.2588, 0.9659}, {-1., 0., 0.}, {0., -0.9659, 0.2588}}}); // rotation matrix from body-fixed frame to sensor-fixed frame
    auto drift_rate = [](double t){ return t * Vector3d({1.0E-8, 1.0E-8, 1.0E-8}); }; // in the sensor-fixed frame
    std::shared_ptr<Gyroscope> gyro = std::make_shared<Gyroscope>(rot_bff2sff, cov_mat_gyro, drift_rate);
    sensors.push_back(std::move(gyro));


    double drag_coeff = 2.1;

    // Create a spacecraft object. The name is used to save the output of the simulation in a file "name.dat".
    Spacecraft sc("sc", initial_state, mass, drag_coeff, specular_reflection_coeff, diffuse_reflection_coeff, inertia_matrix, face_normals, face_areas, cop_positions, residual_dipole, sensors);

    // Integrate the equations of motion.
    sim.integrate(sc, env_model, epoch, t0, t1, dt);
    std::cout << "The results have been saved in \"" << sc.name() << ".dat\"." << std::endl;

    return 0;
}