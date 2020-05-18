//
// Created by Sylvain  on 3/31/20.
//

#ifndef SATOPS_SIM_H
#define SATOPS_SIM_H

#include <string>
#include "dimensions.h"

class Spacecraft;
class EnvironmentModel;

/** Base class used to load the SPICE kernels and run the simulation. */
class Sim {
public:
    explicit Sim(std::string meta_kernel);
    ~Sim();

    /** Integrates the equations of motion. Currently uses a Runge-Kutta 4 integration scheme with fixed time step.
     *
     * @param sc            Reference to the spacecraft
     * @param env_model     Environment model
     * @param epoch         Initial epoch with format "YYYY-MM-DD hh:mm:ss"
     * @param t0            Initial time (sec)
     * @param t1            Final time (sec)
     * @param dt            Integration step (sec)
     */
    void integrate(Spacecraft& sc, const EnvironmentModel& env_model, const std::string& epoch, const dimension::Time t0, const dimension::Time t1, const dimension::Time dt);

private:
    std::string m_meta_kernel;
};


#endif //SATOPS_SIM_H
