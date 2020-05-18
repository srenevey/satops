//
// Created by Sylvain  on 4/1/20.
//

#ifndef SATOPS_VELOCITY_H
#define SATOPS_VELOCITY_H

#include "BaseDimension.h"

namespace dimension {
    /** Definition of velocity quantities. The base unit used throughout the program is km/s. */
    class Velocity: public BaseDimension {
    public:
        /**
         * @param velocity in kilometers per second
         */
        explicit Velocity(double velocity = 0.);
        explicit Velocity(const BaseDimension& a);
    };
}

namespace unit {
    // Base unit for velocities used throughout the program is km/s.
    using namespace dimension;
    /** Creates a velocity in milimeters per second (mm/s). */
    Velocity operator "" _mms(long double d);
    /** Creates a velocity in centimeters per second (cm/s). */
    Velocity operator "" _cms(long double d);
    /** Creates a velocity in meters per second (m/s). */
    Velocity operator "" _ms(long double d);
    /** Creates a velocity in kilometers per second (km/s). */
    Velocity operator "" _kms(long double d);
    /** Creates a velocity in kilometers per hour (km/h). */
    Velocity operator "" _kph(long double d);
    /** Creates a velocity in miles per second (mi/s). */
    Velocity operator "" _mps(long double d);
    /** Creates a velocity in miles per hour (mi/h). */
    Velocity operator "" _mph(long double d);
}

#endif //SATOPS_VELOCITY_H
