//
// Created by Sylvain  on 4/6/20.
//

#include "Area.h"
#include "Distance.h"
#include "SatOps/constants.h"

namespace dimension {
    Area::Area(double area): BaseDimension(area) {}
    Area::Area(const BaseDimension& a): BaseDimension(a) {}

    Area operator*(const Distance& d1, const Distance& d2) {
        return Area(d1.data() * d2.data());
    }
}

namespace unit {
    // Base unit for areas used throughout the program is km^2.
    using namespace dimension;
    Area operator "" _cm2(long double a) {
        return Area(a * constants::CM_TO_KM * constants::CM_TO_KM);
    }

    Area operator "" _m2(long double a) {
        return Area(a * constants::M_TO_KM * constants::M_TO_KM);
    }

    Area operator "" _km2(long double a) {
        return Area(a);
    }

    Area operator "" _in2(long double a) {
        return Area(a * constants::IN_TO_KM * constants::IN_TO_KM);
    }

    Area operator "" _ft2(long double a) {
        return Area(a * constants::FT_TO_KM * constants::FT_TO_KM);
    }

    Area operator "" _acre(long double a) {
        return Area(a * 1./640 * constants::MI_TO_KM * constants::MI_TO_KM);
    }

    Area operator "" _ha(long double a) {
        return Area(a * 0.01);
    }

    Area operator "" _mi2(long double a) {
        return Area(a * constants::MI_TO_KM * constants::MI_TO_KM);
    }
}