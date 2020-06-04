//
// Created by Sylvain  on 5/20/20.
//

#include "gtest/gtest.h"
#include "SatOps/Body.h"
#include "SatOps/constants.h"
#include "SatOps/EnvironmentModel.h"
#include "SatOps/ReferenceFrame.h"
#include "SatOps/Sim.h"
#include "SatOps/StateVector.h"
#include "SatOps/Vector.h"
#include <filesystem>
#include <string>

extern "C" {
#include "SpiceUsr.h"
};

using namespace unit;

// The expected values have been retrived from JPL's Horizon system (https://ssd.jpl.nasa.gov/horizons.cgi).
TEST(TestEnvModel, BodyPosition1) {
    std::filesystem::path meta_kernel("/Users/sylvain/Developer/satops/assets/kernels.tm");
    Sim sim(meta_kernel);

    EnvironmentModel env_model(Body::Earth);
    std::string epoch("2020-05-20 15:20:59.9997 TDB");
    double et;
    str2et_c(epoch.c_str(), &et);
    Vector3<dimension::Distance> pos = env_model.getBodyPosition(Body::Sun, et, ReferenceFrame::J2000);
    Vector3<dimension::Distance> exp({7.619191141417447E+07, 1.200460039540816E+08, 5.203913474471303E+07});
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(pos[i], exp[i], 1E-3);
    }
}

TEST(TestEnvModel, BodyPosition2) {
    std::filesystem::path meta_kernel("/Users/sylvain/Developer/satops/assets/kernels.tm");
    Sim sim(meta_kernel);

    EnvironmentModel env_model(Body::MarsBarycenter);
    std::string epoch("2018-12-24 23:34:11.9997 TDB");
    double et;
    str2et_c(epoch.c_str(), &et);
    Vector3<dimension::Distance> pos = env_model.getBodyPosition(Body::JupiterBarycenter, et, ReferenceFrame::ECLIPJ2000);
    Vector3<dimension::Distance> exp({-4.983014599811100E+08, -8.621936597024153E+08, 1.180828806358087E+07});
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(pos[i], exp[i], 1E-3);
    }
}