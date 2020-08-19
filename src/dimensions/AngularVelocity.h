//
// Created by Sylvain  on 4/1/20.
//

#ifndef SATOPS_ANGULARVELOCITY_H
#define SATOPS_ANGULARVELOCITY_H

#include "BaseDimension.h"

namespace dimension {
    /** Definition of angular velocity quantities. The base unit used throughout the program is rad/s. */
    class AngularVelocity : public BaseDimension {
    public:
        /**
         * @param angular_velocity in radians per second
         */
        explicit AngularVelocity(double angular_velocity = 0.);
        explicit AngularVelocity(const BaseDimension& a);
    };
}

namespace unit {
    // Base unit for angular velocities used throughout the program is rad/s.
    using namespace dimension;

    /** Creates an angular velocity in radians per second (rad/s). */
    AngularVelocity operator "" _radps(long double d);
    /** Creates an angular velocity in degrees per second (deg/s). */
    AngularVelocity operator "" _degps(long double d);
}


#endif //SATOPS_ANGULARVELOCITY_H
