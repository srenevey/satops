//
// Created by Sylvain  on 4/19/20.
//

#ifndef SATOPS_MASS_H
#define SATOPS_MASS_H

#include "BaseDimension.h"

namespace dimension {
    /** Definition of mass quantities. The base unit used throughout the program is kg. */
    class Mass : public BaseDimension {
    public:
        explicit Mass(double mass = 0.);
        explicit Mass(const BaseDimension& a);
    };
}

namespace unit {
    // Base unit for masses used throughout the program is kg.
    using namespace dimension;
    /** Creates a mass in grams (g). */
    Mass operator "" _g(long double m);
    /** Creates a mass in kilograms (kg). */
    Mass operator "" _kg(long double m);
    /** Creates a mass in metric tons. */
    Mass operator "" _t(long double m);
    /** Creates a mass in ounces (oz). */
    Mass operator "" _oz(long double m);
    /** Creates a mass in pounds (lbs). */
    Mass operator "" _lb(long double m);
    /** Creates a mass in US tons. */
    Mass operator "" _ust(long double m);
}


#endif //SATOPS_MASS_H
