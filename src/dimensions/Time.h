//
// Created by Sylvain  on 4/1/20.
//

#ifndef SATOPS_TIME_H
#define SATOPS_TIME_H

#include "BaseDimension.h"

namespace dimension {
    /** Definition of time quantities. The base unit used throughout the program is second. */
    class Time: public BaseDimension {
    public:
        /**
         * @param time in seconds
         */
        explicit Time(double time = 0.);
        explicit Time(const BaseDimension& a);
    };
}

namespace unit {
    // Base unit for time used throughout the program is second.
    using namespace dimension;
    /** Creates a time in seconds (s). */
    Time operator "" _s(unsigned long long t);
    /** Creates a time in seconds (s). */
    Time operator "" _s(long double t);
    /** Creates a time in minutes (min). */
    Time operator "" _min(unsigned long long t);
    /** Creates a time in minutes (min). */
    Time operator "" _min(long double t);
    /** Creates a time in hours (h). */
    Time operator "" _h(unsigned long long t);
    /** Creates a time in hours (h). */
    Time operator "" _h(long double t);
    /** Creates a time in days (d). */
    Time operator "" _d(unsigned long long t);
    /** Creates a time in days (d). */
    Time operator "" _d(long double t);
}

#endif //SATOPS_TIME_H
