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

TEST(TestStateVector, J2000TransformToFrameITRF93) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);
    StateVector state = createState();
    StateVector s_itrf93 = state.transformToFrame(ReferenceFrame::ITRF93);

    // Expected
    std::string epoch = "2020-04-05 16:00:00";
    Vector3<dimension::Distance> position(ReferenceFrame::ITRF93, {-1569.892782880579_km, 6836.720481416324_km, -1891.443018014624_km});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::ITRF93, {-6.240368200367_kms, -1.791615423803_kms, 3.265717969670_kms});
    Vector3<dimension::AngularVelocity> ang_velocity(ReferenceFrame::BODY, {0., 0., 0.});
    StateVector exp_itrf93(epoch, position, velocity, Quaternion::identity(), ang_velocity);
    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_NEAR(s_itrf93[i], exp_itrf93[i], 1E-12);
    }

}

TEST(TestStateVector, ITRF93TransformToFrameJ2000) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);
    std::string epoch("2018-12-21 09:36:55");
    Vector3<dimension::Distance> position(ReferenceFrame::ITRF93, {7394.13194_km, -9837.143_km, 12802.1259_km});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::ITRF93, {5.9230_kms, 3.519238_kms, -6.098234_kms});
    Vector3<dimension::AngularVelocity> ang_velocity(ReferenceFrame::BODY, {0., 0., 0.});
    StateVector state(epoch, position, velocity, Quaternion::identity(), ang_velocity);
    StateVector state_j2000 = state.transformToFrame(ReferenceFrame::J2000);

    //Expected
    position = Vector3<dimension::Distance>(ReferenceFrame::J2000, {-12281.663765472878367_km, -178.580152107684768_km, 12824.423627951931849_km});
    velocity = Vector3<dimension::Velocity>(ReferenceFrame::J2000, {-0.644454159455230_kms, -7.756353603465013_kms, -6.097292618019948_kms});
    StateVector exp_j2000(epoch, position, velocity, Quaternion::identity(), ang_velocity);
    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_NEAR(state_j2000[i], exp_j2000[i], 1E-12);
    }
}

TEST(TestStateVector, ECLIPJ2000TransformToFrameJ2000) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);
    std::string epoch = "2012-02-28 11:21:59";
    Vector3<dimension::Distance> position(ReferenceFrame::ECLIPJ2000, {5921.132938_km, 22850.33480_km, -11803.4363});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::ECLIPJ2000, {-2.59108_kms, 7.2841841_kms, -3.523823_kms});
    Vector3<dimension::AngularVelocity> ang_velocity(ReferenceFrame::BODY, {0., 0., 0.});
    StateVector s_eclipj2000(epoch, position, velocity, Quaternion::identity(), ang_velocity);
    StateVector s_j2000 = s_eclipj2000.transformToFrame(ReferenceFrame::J2000);

    // Expected
    position = Vector3<dimension::Distance>(ReferenceFrame::J2000, {5921.132937999999740_km, 25659.909612912695593_km, -1740.099887190197933_km});
    velocity = Vector3<dimension::Velocity>(ReferenceFrame::J2000, {-2.591080000000000_kms, 8.084804539507012_kms, -0.335562357824344_kms});
    ang_velocity = Vector3<dimension::AngularVelocity>(ReferenceFrame::BODY, {0., 0., 0.});
    StateVector exp1(epoch, position, velocity, Quaternion::identity(), ang_velocity);
    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_NEAR(s_j2000[i], exp1[i], 1E-12);
    }
}

TEST(TestStateVector, J2000TrasformToFrameECLIPJ2000) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);
    std::string epoch = "1995-10-05 21:36:15";
    Vector3<dimension::Distance> position(ReferenceFrame::J2000, {-4657.3245_km, -9421.13421_km, 36421.49384_km});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::J2000, {7.942381_kms, -6.194518_kms, 8.23485319_kms});
    Vector3<dimension::AngularVelocity> ang_velocity(ReferenceFrame::BODY, {0., 0., 0.});
    StateVector s_j2000(epoch, position, velocity, Quaternion::identity(), ang_velocity);
    StateVector s_eclipj2000 = s_j2000.transformToFrame(ReferenceFrame::ECLIPJ2000);

    // Expected
    position = Vector3<dimension::Distance>(ReferenceFrame::J2000, {-4657.324499999999716_km, 5843.916592445603783_km, 37163.579243669861171_km});
    velocity = Vector3<dimension::Velocity>(ReferenceFrame::J2000, {7.942381000000000_kms, -2.407722666729617_kms, 10.019367838007227_kms});
    ang_velocity = Vector3<dimension::AngularVelocity>(ReferenceFrame::BODY, {0., 0., 0.});
    StateVector exp_eclipj2000(epoch, position, velocity, Quaternion::identity(), ang_velocity);
    for (std::size_t i = 0; i < 6; ++i) {
        EXPECT_NEAR(s_eclipj2000[i], exp_eclipj2000[i], 1E-12);
    }
}

TEST(TestStateVector, NormInf) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);
    StateVector state = createState();
    EXPECT_DOUBLE_EQ(state.normInf(), 7009.22977);
}

TEST(TestStateVector, NormalizeQuaternion) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);
    StateVector state = createState();
    state.normalizeQuaternion();
    Quaternion q = state.orientation();
    Quaternion exp(ReferenceFrame::BODY, 0.267600421808970, 0.700483457088185, 0.354177028864813, 0.558812645542260);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(q[i], exp[i], 1E-12);
    }
}

TEST(TestStateVector, Abs) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);

    std::string epoch = "2018-06-30 17:08:00";
    Vector3<dimension::Distance> position(ReferenceFrame::J2000, {-7009.22977_km, 356.45_km, -1877.9_km});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::J2000, {0.0_kms, -7.0033_kms, 3.2657_kms});
    Quaternion orientation(ReferenceFrame::J2000, 0.34, -0.89, 0.45, 0.71);
    Vector3<dimension::AngularVelocity> angular_velocity(ReferenceFrame::BODY, {-0.2_degs, 0.0_degs, -0.1_degs});
    StateVector state(epoch, position, velocity, orientation, angular_velocity);
    StateVector state_abs = abs(state);

    position = Vector3<dimension::Distance>(ReferenceFrame::J2000, {7009.22977_km, 356.45_km, 1877.9_km});
    velocity = Vector3<dimension::Velocity>(ReferenceFrame::J2000, {0.0_kms, 7.0033_kms, 3.2657_kms});
    orientation = Quaternion(ReferenceFrame::J2000, 0.34, 0.89, 0.45, 0.71);
    angular_velocity = Vector3<dimension::AngularVelocity>(ReferenceFrame::BODY, {0.2_degs, 0.0_degs, 0.1_degs});
    StateVector exp(epoch, position, velocity, orientation, angular_velocity);
    for (std::size_t i = 0; i < 13; ++i) {
        EXPECT_DOUBLE_EQ(state_abs[i], exp[i]);
    }
}

TEST(TestStateVector, ComputeGeodeticCoord) {
    std::string meta_kernel(R"(/Users/sylvain/Developer/satops/assets/kernels.tm)");
    Sim sim(meta_kernel);

    std::string epoch = "2020-05-21 16:18:24";
    Vector3<dimension::Distance> position(ReferenceFrame::ITRF93, {-7009.22977_km, 356.45_km, -1877.9_km});
    Vector3<dimension::Velocity> velocity(ReferenceFrame::ITRF93, {0.0_kms, -7.0033_kms, 3.2657_kms});
    Quaternion orientation(ReferenceFrame::J2000, 0.34, -0.89, 0.45, 0.71);
    Vector3<dimension::AngularVelocity> angular_velocity(ReferenceFrame::BODY, {-0.2_degs, 0.0_degs, -0.1_degs});
    StateVector state(epoch, position, velocity, orientation, angular_velocity);
    Vector3d geodetic_coord = state.computeGeodeticCoord();

    Vector3d exp({-0.262923108342161, 3.090782049883596, 8.884796039663423E+02});
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(geodetic_coord[i], exp[i], 1E-12);
    }

    position = Vector3<dimension::Distance>(ReferenceFrame::ITRF93, {13921.21449_km, -2481.4194_km, -11853.3492_km});
    state.setPosition(position);
    geodetic_coord = state.computeGeodeticCoord();
    exp = Vector3d({-0.698772490083819, -0.176394761342321, 12082.237613403966495});
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(geodetic_coord[i], exp[i], 1E-12);
    }
}