//
// Created by Sylvain  on 5/14/20.
//

#include "gtest/gtest.h"
#include "SatOps/StateVector.h"
#include "SatOps/ReferenceFrame.h"
#include "SatOps/dimensions.h"
#include "SatOps/Quaternion.h"
#include "SatOps/Vector.h"
#include "SatOps/Sim.h"
#include <cstddef>
#include <string>


using namespace unit;

StateVector createState() {
    std::string epoch = "2020-04-05 16:00:00";
    Vector3<dimension::Distance> position(ReferenceFrame::J2000, {-7009.22977_km, 356.45_km, -1877.9_km});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::J2000, {0.0_kms, -7.0033_kms, 3.2657_kms});
    Quaternion orientation(ReferenceFrame::J2000, 0.34, 0.89, 0.45, 0.71);
    Vector3<dimension::AngularVelocity> angular_velocity(ReferenceFrame::BODY, {0.2_degs, 0.0_degs, 0.1_degs});
    return StateVector(epoch, position, velocity, orientation, angular_velocity);
}

TEST(TestStateVector, TransformToFrame) { // TODO: implement

}

TEST(TestStateVector, NormInf) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/SatOps/assets/kernels.tm)");
    Sim sim(meta_kernel);
    StateVector state = createState();
    EXPECT_DOUBLE_EQ(state.normInf(), 7009.22977);
}

TEST(TestStateVector, NormalizeQuaternion) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/SatOps/assets/kernels.tm)");
    Sim sim(meta_kernel);
    StateVector state = createState();
    state.normalizeQuaternion();
    Quaternion q = state.orientation();
    Quaternion exp(ReferenceFrame::BODY, 0.267600421808970, 0.700483457088185, 0.354177028864813, 0.558812645542260);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(q[i], exp[i], 1E-12);
    }
}