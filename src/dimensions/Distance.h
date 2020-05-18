//
// Created by Sylvain  on 4/1/20.
//

#ifndef SATOPS_DISTANCE_H
#define SATOPS_DISTANCE_H

#include "BaseDimension.h"
#include "SatOps/constants.h"

namespace dimension {
    class Area;

    /** Definition of distance quantities. The base unit used throughout the program is km. */
    class Distance: public BaseDimension {
    public:
        /**
         * @param distance in kilometers
         */
        explicit Distance(double distance = 0.);
        explicit Distance(const BaseDimension& a);
        friend Area operator*(const Distance& d1, const Distance& d2);
    };
}

namespace unit {
    // Base unit for distances used throughout the program is km.
    using namespace dimension;
    /** Creates a distance in milimeters (mm). */
    Distance operator "" _mm(long double d);
    /** Creates a distance in centimeters (cm). */
    Distance operator "" _cm(long double d);
    /** Creates a distance in decimeters (dm). */
    Distance operator "" _dm(long double d);
    /** Creates a distance in meters (m). */
    Distance operator "" _m(long double d);
    /** Creates a distance in kilometers (km). */
    Distance operator "" _km(long double d);
    /** Creates a distance in astronomical units (AU). */
    Distance operator "" _au(long double d);
    /** Creates a distance in inches (in). */
    Distance operator "" _in(long double d);
    /** Creates a distance in feet (ft). */
    Distance operator "" _ft(long double d);
    /** Creates a distance in mile (mi)s. */
    Distance operator "" _mi(long double d);
}


#endif //SATOPS_DISTANCE_H
