//
// Created by Sylvain  on 3/31/20.
//

#include "Integrator.h"
#include "SatOps/Spacecraft.h"
#include "EqMotion.h"
#include <iostream>

void Integrator::integrate(Spacecraft &sc, const EqMotion& eom, const double t0, const double t1, const double dt) {
    //boost::numeric::odeint::integrate_adaptive( boost::numeric::odeint::make_controlled<stepper>(1E-10, 1E-10), eom , sc.m_state, t0, t1, dt, sc);
    std::cout << "Simulation running..." << std::endl;
    boost::numeric::odeint::integrate_adaptive(stepper(), eom , sc.m_state, t0, t1, dt, sc);
}