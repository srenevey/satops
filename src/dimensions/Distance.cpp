//
// Created by Sylvain  on 4/1/20.
//

#include "Distance.h"

namespace dimension {
    Distance::Distance(double distance): BaseDimension(distance) {}
    Distance::Distance(const BaseDimension& a): BaseDimension(a) {}
}

namespace unit {
    // Base unit for distances used throughout the program is km.
    using namespace dimension;
    Distance operator "" _mm(long double d) {
        return Distance(constants::MM_TO_KM * d);
    }

    Distance operator "" _cm(long double d) {
        return Distance(constants::CM_TO_KM * d);
    }

    Distance operator "" _dm(long double d) {
        return Distance(constants::DM_TO_KM * d);
    }

    Distance operator "" _m(long double d) {
        return Distance(constants::M_TO_KM * d);
    }

    Distance operator "" _km(long double d) {
        return Distance(d);
    }

    Distance operator "" _au(long double d) {
        return Distance(constants::AU_TO_KM * d);
    }

    Distance operator "" _in(long double d) {
        return Distance(constants::IN_TO_KM * d);
    }

    Distance operator "" _ft(long double d) {
        return Distance(constants::FT_TO_KM * d);
    }

    Distance operator "" _mi(long double d) {
        return Distance(constants::MI_TO_KM * d);
    }
}