//
// Created by Sylvain  on 4/19/20.
//

#include "Mass.h"
#include "SatOps/constants.h"

namespace dimension {
    Mass::Mass(double mass): BaseDimension(mass) {}
    Mass::Mass(const BaseDimension& a): BaseDimension(a) {}
}

namespace unit {
    using namespace dimension;

    Mass operator "" _g(long double m) {
        return Mass(constants::G_TO_KG * m);
    }

    Mass operator "" _kg(long double m) {
        return Mass(m);
    }

    Mass operator "" _t(long double m) {
        return Mass(constants::T_TO_KG * m);
    }

    Mass operator "" _oz(long double m) {
        return Mass(constants::OZ_TO_KG * m);
    }

    Mass operator "" _lb(long double m) {
        return Mass(constants::LB_TO_KG * m);
    }

    Mass operator "" _ust(long double m) {
        return Mass(constants::UST_TO_KG * m);
    }
}