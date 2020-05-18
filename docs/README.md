## SatOps

SatOps is a library that propagates the trajectory and attitude of a spacecraft subject to various orbital perturbations. Currently supported are atmospheric drag (exponentional model or EarthGRAM 2016), geopotential perturbations, third-body effects, solar radiation pressure, and magnetic perturbations. An example is provided in *examples/cubesat.cpp*.

Note: the propagator is still a work in progress and has not been validated yet. The output might be inaccurate / incorrect.

## Prerequisites

The C++ propagator uses the EGM2008 geopotential model, the EarthGRAM2016 atmospheric model, the IGRF magnetic model, as well as several kernels used with the Naif SPICE toolbox for the planetary ephemerides. These models should be downloaded before compiling the library.

### Geopotential model
The EGM2008 model can be downloaded [here](http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/first_release.html). Download the model *Spherical Harmonic Coefficients for Earth's Gravitational Potential - "Tide Free" system* and then extract it. Place the file *EGM2008_to2190_TideFree* into *assets*.

```
SatOps
└── assets
    |   EGM2008_to2190_TideFree
|   docs
|   examples
|   extern
|   include
|   LICENSE
|   README.md
|   src
|   tests
```

### Atmosphere model
EarthGRAM2016 can be requested from the [NASA Software Catalog](https://software.nasa.gov/software/MFS-32780-2). Unzip the archive and place it in the *extern* folder

```
SatOps
└── assets
    |   EGM2008_to2190_TideFree 	
|   docs
|   examples
└── extern
    |   earthGRAM2016
|   include
|   LICENSE
|   README.md
|   src
|   tests
```

Copy *NameRef_Linux.txt* from *data/IOfiles* into *earthGRAM2016* and rename it to *NameRef.txt*. Open the file and modifiy as follows:

    atmpath = <ABSOLUTE_PATH_TO_PROJECT_ROOT>/extern/earthGRAM2016/data/IOfiles/
    NCEPpath = <ABSOLUTE_PATH_TO_PROJECT_ROOT>/extern/earthGRAM2016/data/NCEPdata/FixedBin/
    trapath = null
    prtpath = null
    nprpath =  null
    conpath =  null
    rrapath = <ABSOLUTE_PATH_TO_PROJECT_ROOT>/extern/earthGRAM2016/data/RRAdata/
    ...
    iurra = 0
    ...
    ibltest = 0

where <ABSOLUTE_PATH_TO_PROJECT_ROOT> must be modified accordingly. By default, EarthGRAM2016 prompts the user to enter the path to *NameRef.txt* which is not ideal for a seamless integration. The file *earthgram.patch* modifies some of the source files to suppress any prompt and output. To apply the patch, copy the file *earthgram.patch* from *assets* into *extern/earthGRAM2016/src* and then in a terminal, navigate to *extern/earthGRAM2016/src* and type

    patch -p1 < earthgram.patch

The next step is to compile EarthGRAM2016 into a library. Assuming that [CMake](https://cmake.org/) is installed, create a CMakeLists file in *earthGRAM2016/src* and add the following lines:

```cmake
cmake_minimum_required(VERSION 3.14)
project(earthGRAM2016)

set(CMAKE_CXX_STANDARD 20)

add_library(earthGRAM2016 STATIC
    Atmod1.cpp
    Atmod1.h
    AuxProf.cpp
    AuxProf.h
    HWM.cpp
    HWM.h
    Init.cpp
    Init.h
    InitP.cpp
    InitP.h
    JB2008.cpp
    JB2008.h
    Map.cpp
    Map.h
    MET.cpp
    MET.h
    MSIS.cpp
    MSIS.h
    NCEP.cpp
    NCEP.h
    Pert.cpp
    Pert.h
    RRA.cpp
    RRA.h)

set_property(TARGET earthGRAM2016 PROPERTY POSITION_INDEPENDENT_CODE ON)
```

Open a terminal and navigate to *earthGRAM2016/src*. Then, type

```zsh
mkdir build
cd build
cmake ..
cmake --build .
```

The library *earthGRAM2016.a* will be crated in the *build* directory.

### Planetary ephemerides
The propagator uses JPL's ephemerides to retrieve the celestial bodies positions. This process is facilitated by the SPICE toolbox developed by NASA's Naif. The kernels used by the SPICE toolbox can be accessed from the [JPL's Naif website](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/) by going to Data > Generic Kernels > Generic Kernels and navigating in the subfolders. The required kernels are (right click > "Save Link As..." to download the file directly):
- [de430.bsp](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp)
- [earth_latest_high_prec.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc)
- [gm_de431.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc)
- [naif0012.tls](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls)
- [pck00010.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc)

These kernels should be placed in a *kernels* folder inside the *assets* directory.

```
SatOps
└── assets
    |   EGM2008_to2190_TideFree 
    └── kernels
        |   de430.bsp
        |   earth_latest_high_prec.bpc
        |   gm_de431.tpc
        |   naif0012.tls
        |   pck00010.tpc
|   docs
|   examples
└── extern
    |   earthGRAM2016
|   include
|   LICENSE
|   README.md
|   src
|   tests
```

These kernels will be loaded by SPICE through a meta-kernel. The following content should be put in a file named *kernels.tm* and placed in the *assets* folder. Modify <ABSOLUTE_PATH_TO_PROJECT_ROOT> to reflect your project's root directory path:

    KPL/MK
       \begindata
       PATH_VALUES     = ( '<ABSOLUTE_PATH_TO_PROJECT_ROOT>' )
       
       PATH_SYMBOLS    = ( 'ROOT' )
    
       KERNELS_TO_LOAD = (  '$ROOT/assets/kernels/de430.bsp',
                            '$ROOT/assets/kernels/earth_latest_high_prec.bpc',
                            '$ROOT/assets/kernels/gm_de431.tpc',
                            '$ROOT/assets/kernels/naif0012.tls',
                            '$ROOT/assets/kernels/pck00010.tpc')
     
       \begintext

### Geomagnetic Model

The International Geomagnetic Reference Field (IGRF) is used to model the Earth magnetic field. The text file containing the coefficients can be downloaded [here](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) and should be placed in the *assets* folder alongside the other models. The IGRF-13 version is provided in this repository.


## Compilation

The SatOps library is based on boost and CSPICE. These libraries should be downloaded and installed before compiling the library. 

Create a build directory to keep the project directory clean:

```zsh
mkdir build
cd build
cmake ..
```

To compile the library, the examples, and the documentation (requires doxygen):

```zsh
cmake --build . 
```

To compile the library only:

```zsh
cmake --build . --target SatOps
```

To compile a specific example:

```zsh
cmake --build . --target cubesat
```

To compile the documentation:

```zsh
cmake --build . --target doc
```

To compile the tests (googletest must be installed in *extern* prior to building the tests):

```zsh
mkdir build
cd build
cmake .. -DBUILD_TESTS=ON
cmake --build .
```

