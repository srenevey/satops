//
// Created by Sylvain  on 4/1/20.
//

#include "Time.h"
#include "SatOps/constants.h"

namespace dimension {
    Time::Time(const double time): BaseDimension(time) {}
    Time::Time(const BaseDimension& a): BaseDimension(a) {}
}

namespace unit {
    using namespace dimension;
    Time operator "" _s(unsigned long long t) {
        return Time(t);
    }

    Time operator "" _s(long double t) {
        return Time(t);
    }

    Time operator "" _min(unsigned long long t) {
        return Time(60.0 * t);
    }

    Time operator "" _min(long double t) {
        return Time(60.0 * t);
    }

    Time operator "" _h(unsigned long long t) {
        return Time(t * constants::HOURS_TO_SEC);
    }

    Time operator "" _h(long double t) {
        return Time(t * constants::HOURS_TO_SEC);
    }

    Time operator "" _d(unsigned long long t) {
        return Time(24 * constants::HOURS_TO_SEC * t);
    }

    Time operator "" _d(long double t) {
        return Time(24 * constants::HOURS_TO_SEC * t);
    }
}