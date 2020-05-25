//
// Created by Sylvain  on 5/13/20.
//

#include "gtest/gtest.h"
#include <cstddef>
#include "SatOps/Matrix.h"

TEST(TestMatrix, Size) {
    Matrix<double, 4, 3> mat({{{3.1, 4., -9.2}, {-3.2, 1.65, 3.2}, {0.8, -7.5, 3.7}, {-10.4, 6.42, 1.4}}});
    std::pair<std::size_t, std::size_t> dim = mat.size();
    EXPECT_EQ(dim.first, 4);
    EXPECT_EQ(dim.second, 3);
}

TEST(TestMatrix, ScalarMul) {
    Matrix<double, 4, 3> m1({{{3.1, 4., -9.2}, {-3.2, 1.65, 3.2}, {0.8, -7.5, 3.7}, {-10.4, 6.42, 1.4}}});
    Matrix<double, 4, 3> prod1 = 4.3 * m1;
    Matrix<double, 4, 3> exp1({{{13.33, 17.2, -39.56}, {-13.76, 7.095, 13.76}, {3.44, -32.25, 15.91}, {-44.72, 27.606, 6.02}}});
    std::pair<std::size_t, std::size_t> dim = exp1.size();
    for (std::size_t i = 0; i < dim.first; ++i) {
        for (std::size_t j = 0; j < dim.second; ++j) {
            EXPECT_DOUBLE_EQ(prod1[i][j], exp1[i][j]);
        }
    }
}

TEST(TestMatrix, Inverse) {
    Matrix<double, 4, 4> m1({{{3.1, 4., -9.2, 14.2}, {-3.2, 1.65, 3.2, -11.8}, {0.8, -7.5, 3.7, 8.9}, {-10.4, 6.42, 1.4, 3.4}}});
    Matrix<double, 4, 4> m1_inv = m1.inverse();
    Matrix<double, 4, 4> exp1({{{-0.636546283809712, -1.186031760389143, -0.570313708686156, 0.035168960239060},
                                {-0.735261871087284, -1.369782130418417, -0.699724594673285, 0.148481859733628},
                                {-0.920338179346337, -1.591246363662943, -0.702924091616203, 0.161211609082193},
                                {-0.179772555339934, -0.386171800480216, -0.133804983609074, 0.054943351259196}}});
    std::pair<std::size_t, std::size_t> dim = exp1.size();
    for (std::size_t i = 0; i < dim.first; ++i) {
        for (std::size_t j = 0; j < dim.second; ++j) {
            EXPECT_NEAR(m1_inv[i][j], exp1[i][j], 1E-12);
        }
    }
}

TEST(TestMatrix, Transpose) {
    Matrix<double, 4, 4> m1({{{3.1, 4., -9.2, 14.2}, {-3.2, 1.65, 3.2, -11.8}, {0.8, -7.5, 3.7, 8.9}, {-10.4, 6.42, 1.4, 3.4}}});
    Matrix<double, 4, 4> m1_transpose = m1.transpose();
    std::pair<std::size_t, std::size_t> dim = m1.size();
    for (std::size_t i = 0; i < dim.first; ++i) {
        for (std::size_t j = 0; j < dim.second; ++j) {
            EXPECT_DOUBLE_EQ(m1[i][j], m1_transpose[j][i]);
        }
    }
}

TEST(TestMatrix, MatMul) {
    Matrix<double, 4, 3> m1({{{3.1, 4., -9.2}, {-3.2, 1.65, 3.2}, {0.8, -7.5, 3.7}, {-10.4, 6.42, 1.4}}});
    Matrix<double, 3, 5> m2({{{6.7, -24.5, 13.7, -8.6, 1.2}, {56.8, 23.5, -59.1, 0.38, 13.8}, {1.54, -24.6, 42.1, -50.9, 44.7}}});
    Matrix<double, 4, 5> prod = m1 * m2;
    Matrix<double, 4, 5> exp1({{{233.802, 244.37, -581.25, 443.14, -352.32},
                                {77.208, 38.455, -6.635, -134.733, 161.97},
                                {-414.942, -286.87, 609.98, -198.06, 62.85},
                                {297.132, 371.23, -462.962, 20.6196, 138.696}}});
    std::pair<std::size_t, std::size_t> dim = prod.size();
    for (std::size_t i = 0; i < dim.first; ++i) {
        for (std::size_t j = 0; j < dim.second; ++j) {
            EXPECT_NEAR(prod[i][j], exp1[i][j], 1E-12);
        }
    }
}

TEST(TestMatrix, Add) {
    Matrix<double, 4, 3> m1({{{3.1, 4., -9.2}, {-3.2, 1.65, 3.2}, {0.8, -7.5, 3.7}, {-10.4, 6.42, 1.4}}});
    Matrix<double, 4, 3> m2({{{6.7, 0.2, -4.1}, {8.5, -12.5, 11.2}, {3.9, -0.5, 1.1}, {-22.53, 7.53, -9.75}}});
    Matrix<double, 4, 3> sum = m1 + m2;
    Matrix<double, 4, 3> exp({{{9.8, 4.2, -13.3}, {5.3, -10.85, 14.4}, {4.7, -8.0, 4.8}, {-32.93, 13.95, -8.35}}});
    std::pair<std::size_t, std::size_t> dim = sum.size();
    for (std::size_t i = 0; i < dim.first; ++i) {
        for (std::size_t j = 0; j < dim.second; ++j) {
            EXPECT_DOUBLE_EQ(sum[i][j], exp[i][j]);
        }
    }
}

TEST(TestMatrix, Sub) {
    Matrix<double, 4, 3> m1({{{3.1, 4., -9.2}, {-3.2, 1.65, 3.2}, {0.8, -7.5, 3.7}, {-10.4, 6.42, 1.4}}});
    Matrix<double, 4, 3> m2({{{6.7, 0.2, -4.1}, {8.5, -12.5, 11.2}, {3.9, -0.5, 1.1}, {-22.53, 7.53, -9.75}}});
    Matrix<double, 4, 3> sub = m1 - m2;
    Matrix<double, 4, 3> exp({{{-3.6, 3.8, -5.1}, {-11.7, 14.15, -8.0}, {-3.1, -7.0, 2.6}, {12.13, -1.11, 11.15}}});
    std::pair<std::size_t, std::size_t> dim = sub.size();
    for (std::size_t i = 0; i < dim.first; ++i) {
        for (std::size_t j = 0; j < dim.second; ++j) {
            EXPECT_DOUBLE_EQ(sub[i][j], exp[i][j]);
        }
    }
}