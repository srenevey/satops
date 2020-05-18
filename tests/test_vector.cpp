//
// Created by Sylvain  on 5/12/20.
//

#include "gtest/gtest.h"
#include "SatOps/Vector.h"
#include "SatOps/Matrix.h"
#include <cstddef>

TEST(TestVector, Add) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<double, 4> v2({1.2, 4.1, -5.3, 0.4});
    Vector<double, 4> sum = v1 + v2;
    Vector<double, 4> exp({4.6, 2.9, 4.5, 12.1});
    for (std::size_t i = 0; i < exp.size(); ++i) {
        EXPECT_DOUBLE_EQ(sum[i], exp[i]);
    }

    Vector<int, 5> v3({7, 0, -3, -9, 1});
    Vector<int, 5> v4({-1, 2, 1, -6, 5});
    Vector<int, 5> sum2 = v3 + v4;
    Vector<int, 5> exp2({6, 2, -2, -15, 6});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(sum2[i], exp2[i]);
    }
}

TEST(TestVector, Sub) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<double, 4> v2({1.2, 4.1, -5.3, 0.4});
    Vector<double, 4> sum = v1 - v2;
    Vector<double, 4> exp({2.2, -5.3, 15.1, 11.3});
    for (std::size_t i = 0; i < exp.size(); ++i) {
        EXPECT_DOUBLE_EQ(sum[i], exp[i]);
    }

    Vector<int, 5> v3({7, 0, -3, -9, 1});
    Vector<int, 5> v4({-1, 2, 1, -6, 5});
    Vector<int, 5> sum2 = v3 - v4;
    Vector<int, 5> exp2({8, -2, -4, -3, -4});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(sum2[i], exp2[i]);
    }
}

TEST(TestVector, ScalarMul) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<double, 4> prod1 = 2.1 * v1;
    Vector<double, 4> exp1({7.14, -2.52, 20.58, 24.57});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_DOUBLE_EQ(prod1[i], exp1[i]);
    }

    Vector<int, 5> v2({7, 0, -3, -9, 1});
    Vector<double, 5> prod2 = -1.8 * v2;
    Vector<double, 5> exp2({-12.6, 0., 5.4, 16.2, -1.8});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_DOUBLE_EQ(prod2[i], exp2[i]);
    }

    Vector<int, 5> v3({7, 0, -3, -9, 1});
    Vector<int, 5> prod3 = 4 * v3;
    Vector<int, 5> exp3({28, 0, -12, -36, 4});
    for (std::size_t i = 0; i < exp3.size(); ++i) {
        EXPECT_EQ(prod3[i], exp3[i]);
    }
}

TEST(TestVector, Mul) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<double, 4> v2({1.2, 4.1, -5.3, 0.4});
    Vector<double, 4> prod1 = v1 * v2;
    Vector<double, 4> exp1({4.08, -4.92, -51.94, 4.68});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_DOUBLE_EQ(prod1[i], exp1[i]);
    }

    Vector<int, 5> v3({7, 0, -3, -9, 1});
    Vector<int, 5> v4({-1, 2, 1, -6, 5});
    Vector<int, 5> prod2 = v3 * v4;
    Vector<int, 5> exp2({-7, 0, -3, 54, 5});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(prod2[i], exp2[i]);
    }
}

TEST(TestVector, ScalarDiv) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<double, 4> div1 = v1 / 2.5;
    Vector<double, 4> exp1({1.36, -0.48, 3.92, 4.68});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_DOUBLE_EQ(div1[i], exp1[i]);
    }

    Vector<int, 5> v2({7, 5, -3, -9, 1});
    Vector<double, 5> div2 = v2 / 7.8;
    Vector<double, 5> exp2({0.897435897435897, 0.641025641025641, -0.384615384615385, -1.153846153846154, 0.128205128205128});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_NEAR(div2[i], exp2[i], 1E-12);
    }
}

TEST(TestVector, Div) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<double, 4> v2({1.2, 4.1, -5.3, 0.4});
    Vector<double, 4> div1 = v1 / v2;
    Vector<double, 4> exp1({2.833333333333333, -0.292682926829268, -1.849056603773585, 29.249999999999996});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_NEAR(div1[i], exp1[i], 1E-12);
    }

    Vector<int, 5> v3({7, 1, -3, -9, 11});
    Vector<int, 5> v4({-2, 2, 1, -6, 5});
    Vector<int, 5> div2 = v3 / v4;
    Vector<int, 5> exp2({-3, 0, -3, 1, 2});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(div2[i], exp2[i]);
    }
}

TEST(TestVector, DotProduct) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<double, 4> v2({1.2, 4.1, -5.3, 0.4});
    double dot1 = v1.dot(v2);
    EXPECT_DOUBLE_EQ(dot1, -48.1);

    Vector<int, 5> v3({7, 1, -3, -9, 11});
    Vector<int, 5> v4({-2, 2, 1, -6, 5});
    double dot2 = v3.dot(v4);
    EXPECT_DOUBLE_EQ(dot2, 94.0);

    Vector<double, 5> v5({3.4, -1.2, 9.8, 11.7, 8.3});
    double dot3 = v3.dot(v5);
    EXPECT_DOUBLE_EQ(dot3, -20.799999999999986);
}

TEST(TestVector, CrossProduct) {
    Vector3d v1({3.4, -1.2, 9.8});
    Vector3d v2({1.2, 4.1, -5.3});
    Vector3d cross1 = v1.cross(v2);
    Vector3d exp1({-33.82, 29.78, 15.379999999999997});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_DOUBLE_EQ(cross1[i], exp1[i]);
    }

    Vector<int, 3> v3({7, 1, -3});
    Vector<int, 3> v4({-2, 2, 1});
    Vector3d cross2 = v3.cross(v4);
    Vector3d exp2({7., -1., 16.});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_DOUBLE_EQ(cross2[i], exp2[i]);
    }

    Vector3d cross3 = v1.cross(v4);
    Vector3d exp3({-20.8, -23., 4.4});
    for (std::size_t i = 0; i < exp3.size(); ++i) {
        EXPECT_DOUBLE_EQ(cross3[i], exp3[i]);
    }
}

TEST(TestVector, Norm) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<int, 3> v2({7, 1, -3});
    EXPECT_DOUBLE_EQ(v1.norm(), 15.682155464093576);
    EXPECT_DOUBLE_EQ(v2.norm(), 7.681145747868608);
}

TEST(TestVector, NormInf) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    Vector<int, 3> v2({7, 1, -3});
    EXPECT_DOUBLE_EQ(v1.normInf(), 11.7);
    EXPECT_DOUBLE_EQ(v2.normInf(), 7.);
}

TEST(TestVector, Normalize) {
    Vector<double, 4> v1({3.4, -1.2, 9.8, 11.7});
    v1.normalize();
    Vector<double, 4> exp1({0.216806931150808, -0.076520093347344, 0.624914095669975, 0.746070910136603});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_NEAR(v1[i], exp1[i], 1E-12);
    }

    Vector<int, 3> v2({7, 1, -3});
    v2.normalize();
    Vector<int, 3> exp2({0, 0, 0});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(v2[i], exp2[i]);
    }
}

TEST(TestVector, Abs) {
    Vector<double, 4> v1({3.4, -1.2, -9.8, 11.7});
    Vector<double, 4> abs1 = v1.abs();
    Vector<double, 4> exp1({3.4, 1.2, 9.8, 11.7});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_DOUBLE_EQ(abs1[i], exp1[i]);
    }

    Vector<int, 5> v2({7, 1, -3, -9, 11});
    Vector<int, 5> abs2 = v2.abs();
    Vector<int, 5> exp2({7, 1, 3, 9, 11});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(abs2[i], exp2[i]);
    }
}

TEST(TestVector, Clamp) {
    Vector<double, 4> v1({3.42, -1.25, -9.81, 11.76});
    v1.clamp(-2., 4.);
    Vector<double, 4> exp1({3.42, -1.25, -2., 4.});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_DOUBLE_EQ(v1[i], exp1[i]);
    }

    Vector<int, 5> v2({7, 1, -3, -9, 11});
    v2.clamp(-4, 8);
    Vector<int, 5> exp2({7, 1, -3, -4, 8});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(v2[i], exp2[i]);
    }
}

TEST(TestVector, RoundDown) {
    Vector<double, 4> v1({3.43, -1.25, -9.81, 11.76});
    v1.roundDown(0.02);
    Vector<double, 4> exp1({3.42, -1.24, -9.80, 11.76});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_DOUBLE_EQ(v1[i], exp1[i]);
    }

    Vector<int, 5> v2({7, 4, -3, -8, 11});
    v2.roundDown(2.);
    Vector<int, 5> exp2({6, 4, -2, -8, 10});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_EQ(v2[i], exp2[i]);
    }
}

TEST(TestVector, Size) {
    Vector<double, 4> v1({3.43, -1.25, -9.81, 11.76});
    Vector<int, 7> v2({5, -2, 5, 0, -1, 10, -94});
    EXPECT_EQ(v1.size(), 4);
    EXPECT_EQ(v2.size(), 7);
}

TEST(TestVector, MatProd) {
    Vector<double, 4> v1({3.43, -1.25, -9.81, 11.76});
    Matrix<double, 3, 4> m1({{{1.2, 3.5, 8.2, -9.1}, {11.3, -21.7, 5.9, 1.1}, {-8.3, 4.2, 0.6, 3.4}}});
    Vector3d prod1 = m1 * v1;
    Vector3d exp1({-187.717, 20.941, 0.379});
    for (std::size_t i = 0; i < exp1.size(); ++i) {
        EXPECT_NEAR(prod1[i], exp1[i], 1E-12);
    }

    Vector<int, 3> v2({5, -2, 3});
    Matrix<double, 4, 3> m2({{{4.5, 1.6, -7.8}, {3.1, -0.8, 2.2}, {-10.3, 5.5, -0.4}, {8.2, 12.2, 1.6}}});
    Vector<double, 4> prod2 = m2 * v2;
    Vector<double, 4> exp2({-4.1, 23.7, -63.7, 21.4});
    for (std::size_t i = 0; i < exp2.size(); ++i) {
        EXPECT_NEAR(prod2[i], exp2[i], 1E-12);
    }
}