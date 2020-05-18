add_executable(cubesat cubesat.cpp)
target_compile_features(cubesat PRIVATE cxx_std_20)
target_link_libraries(cubesat PRIVATE SatOps)