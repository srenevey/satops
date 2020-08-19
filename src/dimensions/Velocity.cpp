//
// Created by Sylvain  on 4/1/20.
//

#include "Velocity.h"
#include "SatOps/constants.h"

namespace dimension {
    Velocity::Velocity(double velocity): BaseDimension(velocity) {}
    Velocity::Velocity(const BaseDimension& a): BaseDimension(a) {}
}

namespace unit {
    // Base unit for velocities used throughout the program is km/s.
    using namespace dimension;
    Velocity operator "" _mmps(long double d) {
        return Velocity(constants::MM_TO_KM * d);
    }

    Velocity operator "" _cmps(long double d) {
        return Velocity(constants::CM_TO_KM * d);
    }

    Velocity operator "" _mps(long double d) {
        return Velocity(constants::M_TO_KM * d);
    }

    Velocity operator "" _kmps(long double d) {
        return Velocity(d);
    }

    Velocity operator "" _kmph(long double d) {
        return Velocity(d / constants::HOURS_TO_SEC);
    }

    Velocity operator "" _mips(long double d) {
        return Velocity(constants::MI_TO_KM * d);
    }

    Velocity operator "" _miph(long double d) {
        return Velocity(constants::MI_TO_KM * d / constants::HOURS_TO_SEC);
    }

}