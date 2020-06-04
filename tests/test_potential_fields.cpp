//
// Created by Sylvain  on 5/30/20.
//

#include "gtest/gtest.h"
#include "SatOps/GravitationalField.h"
#include "SatOps/MagneticField.h"
#include "SatOps/Vector.h"
#include "SatOps/Matrix.h"
#include <filesystem>
#include <iostream>

// Test gravitational and magnetic fields using values provided in Gottlieb, R.G., Fast Gravity, Gravity Partials,
// Normalized Gravity, Gravity Gradient Torque and Magnetic Field: Derivation, Code and Data. NASA Contractor Report 188243, February 1993.
TEST(TestPotentialFields, IGRF13MagField) {
    std::filesystem::path mag_model_path("/Users/sylvain/Developer/satops/assets/igrf13coeffs_modified.txt");
    MagneticField model(MagneticField::IGRF13, 10, mag_model_path, 1985.0);
    Vector3d pos({5489.150, 802.222, 3140.916});
    Matrix3d jacobian;
    Vector3d mag_field = model.getMagneticField(pos, &jacobian);
    Vector3d exp({-3.75259753018348E-05, -6.16002442001094E-06, 1.35121172654619E-05});
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(mag_field[i], exp[i], 1E-12);
    }
}

TEST(TestPotentialFields, GEM10GravityField) {
    std::filesystem::path geopot_model_path("/Users/sylvain/Developer/satops/assets/GEM10");
    GravitationalField model(GravitationalField::GEM10, 4, geopot_model_path);
    Vector3d pos({5489.150, 802.222, 3140.916});
    Matrix3d jacobian;
    Vector3d gravity_field = model.getAcceleration(pos, &jacobian);
    Vector3d exp({-8.44269212018857E-03, -1.233936337854853E-03, -4.84659352346614E-03});
    Matrix3d exp_jacobian({{{1.87779191647821E-6, 4.99270741320439E-7, 1.96515882331155E-6}, {4.99270741320439E-7, -1.46519995493325E-6, 2.87214112813307E-7}, {1.96515882331155E-6, 2.87214112813307E-7, -4.12591961544953E-7}}});

    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(gravity_field[i], exp[i], 1E-12);
    }

    std::pair<std::size_t, std::size_t> dim = jacobian.size();
    for (std::size_t i = 0; i < dim.first; ++i) {
        for (std::size_t j = 0; j < dim.second; ++j) {
            EXPECT_NEAR(jacobian[i][j], exp_jacobian[i][j], 1E-7);
        }
    }
}