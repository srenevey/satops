//
// Created by Sylvain  on 3/31/20.
//

#ifndef SATOPS_INTEGRATOR_H
#define SATOPS_INTEGRATOR_H

#include <boost/numeric/odeint.hpp>

class Spacecraft;
class StateVector;
class EqMotion;

/** Simple wrapper around boost's odeint integration method. */
class Integrator {
public:
    static void integrate(Spacecraft& sc, const EqMotion& eom, double t0, double t1, double dt);

private:
    //typedef boost::numeric::odeint::runge_kutta_dopri5< State, double, State, double, boost::numeric::odeint::vector_space_algebra > stepper;
    typedef boost::numeric::odeint::runge_kutta4< StateVector, double, StateVector, double, boost::numeric::odeint::vector_space_algebra > stepper;
};


#endif //SATOPS_INTEGRATOR_H
