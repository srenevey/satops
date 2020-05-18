//
// Created by Sylvain  on 4/7/20.
//

#ifndef SATOPS_SENSOR_H
#define SATOPS_SENSOR_H

#include "Vector.h"

class StateVector;

/** Base class representing a sensor */
class Sensor {
public:
    Sensor() {}
    virtual ~Sensor() {}

    /**
     * @param state     State vector of the spacecraft.
     * @return          Measurement in the sensor-fixed frame.
     */
    virtual Vector3d measure(const StateVector& state) = 0;
};


#endif //SATOPS_SENSOR_H
