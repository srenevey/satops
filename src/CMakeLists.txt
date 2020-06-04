set(HEADER_LIST
        "${SatOps_SOURCE_DIR}/include/SatOps/Atmosphere.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Body.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/constants.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Dimensions.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/EnvironmentModel.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/GravitationalField.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Gyroscope.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/MagneticField.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Magnetometer.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Matrix.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Quaternion.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/ReferenceFrame.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Sensor.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Sim.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Spacecraft.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/StateVector.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/Vector.h"
        "${SatOps_SOURCE_DIR}/include/SatOps/SatOps.h"
        dimensions/Acceleration.h
        dimensions/AngularAcceleration.h
        dimensions/AngularVelocity.h
        dimensions/Area.h
        dimensions/BaseDimension.h
        dimensions/Distance.h
        dimensions/Mass.h
        dimensions/Time.h
        dimensions/Velocity.h
        EqMotion.h
        Integrator.h
        PotentialField.h
        transformations.h)

set(SRCS_LIST
        dimensions/Acceleration.cpp
        dimensions/AngularAcceleration.cpp
        dimensions/AngularVelocity.cpp
        dimensions/Area.cpp
        dimensions/BaseDimension.cpp
        dimensions/Distance.cpp
        dimensions/Mass.cpp
        dimensions/Time.cpp
        dimensions/Velocity.cpp
        Atmosphere.cpp
        EnvironmentModel.cpp
        EqMotion.cpp
        GravitationalField.cpp
        Gyroscope.cpp
        Integrator.cpp
        MagneticField.cpp
        Magnetometer.cpp
        Quaternion.cpp
        Sim.cpp
        Spacecraft.cpp
        PotentialField.cpp
        StateVector.cpp)

# Make library
add_library(SatOps ${SRCS_LIST} ${HEADER_LIST})
target_include_directories(SatOps PUBLIC ../include ../extern/earthGRAM2016/src /usr/local/cspice/include "${Boost_INCLUDE_DIRS}")
target_link_libraries(SatOps PRIVATE Boost::boost ${EARTHGRAM_LIB} ${CSPICE_LIB} ${CSUPPORT_LIB})
target_compile_features(SatOps PUBLIC cxx_std_20)

# IDEs should put the headers in a nice place
# source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${HEADER_LIST})
