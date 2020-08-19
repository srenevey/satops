//
// Created by Sylvain  on 4/13/20.
//

#include "Acceleration.h"

namespace dimension {
    Acceleration::Acceleration(double acceleration): BaseDimension(acceleration) {}
    Acceleration::Acceleration(const BaseDimension& a): BaseDimension(a) {}
}

namespace unit {
    // Base unit for accelerations used throughout the program is km/s^2.
    using namespace dimension;

    Acceleration operator "" _kmps2(long double a) {
        return Acceleration(a);
    }
}