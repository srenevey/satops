cmake_minimum_required(VERSION 3.11...3.16)

project(SatOps
    VERSION 1.2
    DESCRIPTION "A satellite simulation environment"
    LANGUAGES CXX)

find_package(Boost REQUIRED)
find_library(CSPICE_LIB cspice.a PATHS "/usr/local/cspice/lib")
find_library(CSUPPORT_LIB csupport.a PATHS "/usr/local/cspice/lib")
find_library(EARTHGRAM_LIB libearthGRAM2016.a PATHS "${PROJECT_SOURCE_DIR}/extern/earthGRAM2016/src/build")

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(docs)


if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR SATOPS_BUILD_TESTING) AND BUILD_TESTING)
    add_subdirectory(tests)
endif()


option(BUILD_TESTS "Build the tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
endif()