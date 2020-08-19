//
// Created by Sylvain  on 4/13/20.
//

#ifndef SATOPS_ACCELERATION_H
#define SATOPS_ACCELERATION_H

#include "BaseDimension.h"

namespace dimension {
    /** Definition of acceleration quantities. The base unit used throughout the program is km/s<sup>2</sup>. */
    class Acceleration : public BaseDimension {
    public:
        /**
         * @param acceleration in kilometers per second
         */
        explicit Acceleration(double acceleration = 0.);
        explicit Acceleration(const BaseDimension& a);
    };
}

namespace unit {
    // Base unit for accelerations used throughout the program is km/s^2.
    using namespace dimension;
    /** Creates an acceleration in kilometers per square second (km/s<sup>2</sup>). */
    Acceleration operator "" _kmps2(long double a);
}

#endif //SATOPS_ACCELERATION_H
