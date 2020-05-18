//
// Created by Sylvain  on 5/13/20.
//

#include "gtest/gtest.h"
#include "SatOps/Quaternion.h"
#include "SatOps/Matrix.h"
#include "SatOps/ReferenceFrame.h"

TEST(TestQuaternion, CrossProdMat) {
    Quaternion q1(ReferenceFrame::BODY, 0.2, 0.1, -0.3, 0.4);
    Matrix3d mat = q1.crossProductMatrix();
    Matrix3d exp({{{0., 0.3, 0.1}, {-0.3, 0., -0.2}, {-0.1, 0.2, 0.}}});
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(mat[i][j], exp[i][j]);
        }
    }
}

TEST(TestQuaternion, Xi) {
    Quaternion q1(ReferenceFrame::BODY, 0.2, 0.1, -0.3, 0.4);
    Matrix<double, 4, 3> xi = q1.xi();
    Matrix<double, 4, 3> exp({{{0.4, 0.3, 0.1}, {-0.3, 0.4, -0.2}, {-0.1, 0.2, 0.4}, {-0.2, -0.1, 0.3}}});
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(xi[i][j], exp[i][j]);
        }
    }
}

TEST(TestQuaternion, Psi) {
    Quaternion q1(ReferenceFrame::BODY, 0.2, 0.1, -0.3, 0.4);
    Matrix<double, 4, 3> psi = q1.psi();
    Matrix<double, 4, 3> exp({{{0.4, -0.3, -0.1}, {0.3, 0.4, 0.2}, {0.1, -0.2, 0.4}, {-0.2, -0.1, 0.3}}});
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(psi[i][j], exp[i][j]);
        }
    }
}

TEST(TestQuaternion, AttitudeMatrix) {
    Quaternion q1(ReferenceFrame::BODY, 0.365148371670111, 0.182574185835055, -0.547722557505166, 0.730296743340222);
    Matrix3d A = q1.attitudeMatrix();
    Matrix3d exp({{{0.333333333333333, -0.666666666666667, -0.666666666666667}, {0.933333333333334, 0.133333333333333, 0.333333333333333}, {-0.133333333333333, -0.733333333333334, 0.666666666666667}}});
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_NEAR(A[i][j], exp[i][j], 1E-12);
        }
    }
}

TEST(TestQuaternion, Conj) {
    Quaternion q1(ReferenceFrame::BODY, 0.2, 0.1, -0.3, 0.4);
    Quaternion q1_conj = q1.conjugate();
    Quaternion exp(ReferenceFrame::BODY, -0.2, -0.1, 0.3, 0.4);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(q1_conj[i], exp[i]);
    }
}

TEST(TestQuaternion, Inverse) {
    Quaternion q1(ReferenceFrame::BODY, 0.2, 0.1, -0.3, 0.4);
    Quaternion q1_inv = q1.inverse();
    Quaternion exp(ReferenceFrame::BODY, -0.666666666666667, -0.333333333333333, 1.000000000000000, 1.333333333333334);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(q1_inv[i], exp[i], 1E-12);
    }
}

TEST(TestQuaternion, Identity) {
    Quaternion q1 = Quaternion::identity();
    Quaternion exp(ReferenceFrame::BODY, 0., 0., 0., 1.);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(q1[i], exp[i]);
    }
}

TEST(TestQuaternion, QuaternionProduct) {
    Quaternion q1(ReferenceFrame::BODY, 0.2, 0.1, -0.3, 0.4);
    Quaternion q2(ReferenceFrame::BODY, -0.3, 1.2, 0.8, 0.2);
    Quaternion prod1 = q1 * q2;
    Quaternion prod2 = q1.quaternionProduct(q2);
    Quaternion exp(ReferenceFrame::BODY, -0.52, 0.57, -0.01, 0.26);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(prod1[i], exp[i], 1E-12);
        EXPECT_NEAR(prod2[i], exp[i], 1E-12);
    }
}


