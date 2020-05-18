//
// Created by Sylvain  on 4/14/20.
//

#include "AngularAcceleration.h"

namespace dimension {
    AngularAcceleration::AngularAcceleration(double acceleration): BaseDimension(acceleration) {}
    AngularAcceleration::AngularAcceleration(const BaseDimension& a): BaseDimension(a) {}
}

namespace unit {
    // Base unit for angular accelerations used throughout the program is rad/s^2.
    using namespace dimension;
    AngularAcceleration operator "" _rads2(long double a) {
        return AngularAcceleration(a);
    }
}