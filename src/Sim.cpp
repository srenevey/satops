//
// Created by Sylvain  on 3/31/20.
//

#include "SatOps/Sim.h"
#include "Integrator.h"
#include "EqMotion.h"
#include "SatOps/Spacecraft.h"
extern "C" {
#include "SpiceUsr.h"
};

Sim::Sim(const std::string meta_kernel): m_meta_kernel(meta_kernel) {
    // Load SPICE kernels
    furnsh_c(m_meta_kernel.c_str());
}

Sim::~Sim() {
    // Clean up kernel space
    unload_c(m_meta_kernel.c_str());
}

void Sim::integrate(Spacecraft& sc, const EnvironmentModel& env_model, const std::string& epoch, const dimension::Time t0, const dimension::Time t1, const dimension::Time dt) {
    double et;
    str2et_c(epoch.c_str(), &et);
    EqMotion eom(sc, env_model, et);
    Integrator::integrate(sc, eom, t0.data(), t1.data(), dt.data());
}