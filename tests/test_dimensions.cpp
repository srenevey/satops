//
// Created by Sylvain  on 5/12/20.
//

#include "gtest/gtest.h"
#include "SatOps/dimensions.h"

using namespace unit;

TEST(TestDimensions, Acceleration) {
    EXPECT_DOUBLE_EQ(34.5_kms2, 34.5);
    EXPECT_DOUBLE_EQ(-154.9_kms2, -154.9);
}

TEST(TestDimensions, AngularAcceleration) {
    EXPECT_DOUBLE_EQ(2.31_rads2, 2.31);
    EXPECT_DOUBLE_EQ(-1.65_rads2, -1.65);
}

TEST(TestDimensions, AngularVelocity) {
    EXPECT_DOUBLE_EQ(0.59_rads, 0.59);
    EXPECT_NEAR(30.5_degs, 0.532325421858270, 1E-12);
    EXPECT_NEAR(0.59_rads + 30.5_degs, 1.122325421858271, 1E-12);
}

TEST(TestDimensions, Area) {
    EXPECT_DOUBLE_EQ(125.9_cm2, 125.9E-10);
    EXPECT_DOUBLE_EQ(31.7_m2, 31.7E-6);
    EXPECT_DOUBLE_EQ(2.34_km2, 2.34);
    EXPECT_DOUBLE_EQ(12.89_in2, 8.316112399999999E-9);
    EXPECT_DOUBLE_EQ(163.4_ft2, 1.5180356736E-5);
    EXPECT_DOUBLE_EQ(12.8_acre, 0.051799504712);
    EXPECT_DOUBLE_EQ(120.7_ha, 1.207);
    EXPECT_DOUBLE_EQ(45.2_mi2, 117.06688064912);

    dimension::Area area(45.2_mi2 - 2.34_km2 + 125.9_cm2 + 31.7_m2 - 12.89_in2 - 163.4_ft2 + 12.8_acre - 120.7_ha);
    EXPECT_DOUBLE_EQ(area, 113.5716966777491);
}

TEST(TestDimensions, Distance) {
    EXPECT_DOUBLE_EQ(3.1_mm, 3.1E-6);
    EXPECT_DOUBLE_EQ(27.9_cm, 27.9E-5);
    EXPECT_DOUBLE_EQ(188.2_dm, 188.2E-4);
    EXPECT_DOUBLE_EQ(-485.2_m, -485.2E-3);
    EXPECT_DOUBLE_EQ(89.4_km, 89.4);
    EXPECT_DOUBLE_EQ(0.3_au, 4.487936120999999E+7);
    EXPECT_DOUBLE_EQ(36.78_in, 9.34212E-4);
    EXPECT_DOUBLE_EQ(-214.8_ft, -0.06547104);
    EXPECT_DOUBLE_EQ(22.5_mi, 36.210149999999999);

    dimension::Distance dist(0.3_au - 89.4_km + 22.5_mi - 3.1_mm - 27.9_cm + 188.2_dm - 485.2_m - 36.78_in + 214.8_ft);
    EXPECT_DOUBLE_EQ(dist, 4.487930761802472E+7);
}

TEST(TestDimensions, Mass) {
    EXPECT_DOUBLE_EQ(356.7_g, 0.3567);
    EXPECT_DOUBLE_EQ(44.1_kg, 44.1);
    EXPECT_DOUBLE_EQ(13.6_t, 13.6E+3);
    EXPECT_DOUBLE_EQ(73.8_oz, 2.092194806625);
    EXPECT_DOUBLE_EQ(969.2_lb, 4.39621725004E+2);
    EXPECT_DOUBLE_EQ(29.85_ust, 2.707947225E+4);

    dimension::Mass mass(29.85_ust - 13.6_t + 44.1_kg - 356.7_g + 73.8_oz - 969.2_lb);
    EXPECT_DOUBLE_EQ(mass, 1.308568601980262E+4);
}

TEST(TestDimensions, Time) {
    EXPECT_DOUBLE_EQ(17_s, 17.0);
    EXPECT_DOUBLE_EQ(89.3_s, 89.3);
    EXPECT_DOUBLE_EQ(41_min, 2460.0);
    EXPECT_DOUBLE_EQ(56.71_min, 3402.6);
    EXPECT_DOUBLE_EQ(3_h, 10800.0);
    EXPECT_DOUBLE_EQ(6.715_h, 24174.0);
    EXPECT_DOUBLE_EQ(17_d, 1468800.0);
    EXPECT_DOUBLE_EQ(21.89_d, 1891296.0);

    dimension::Time t(21.89_d - 17_d + 17_s - 89.3_s - 41_min + 56.71_min - 3_h - 6.715_h);
    EXPECT_DOUBLE_EQ(t, 3.883923E+5);
}

TEST(TestDimensions, Velocity) {
    EXPECT_DOUBLE_EQ(14.6_mms, 14.6E-6);
    EXPECT_DOUBLE_EQ(-67.82_cms, -67.82E-5);
    EXPECT_DOUBLE_EQ(188.91_ms, 188.91E-3);
    EXPECT_DOUBLE_EQ(-894.5_kms, -894.5);
    EXPECT_NEAR(80.0_kph, 0.022222222222222, 1E-12);
    EXPECT_DOUBLE_EQ(125.3_mps, 2.01650302E+2);
    EXPECT_NEAR(-89.5_mph, -0.040009980555556, 1E-12);

    dimension::Velocity vel(125.3_mps + 89.5_mph - 894.5_kms + 80.0_kph - 14.6_mms - 67.82_cms + 188.91_ms);
    EXPECT_DOUBLE_EQ(vel, -6.925992485972222E+2);
}