cmake_minimum_required(VERSION 3.4...3.18)
project(rte-rrtmgp LANGUAGES)

# Define flags to identify the operating system
set(LINUX     FALSE)
set(APPLE_ARM FALSE)
set(WINDOWS   FALSE)

# Determine the operating system and architecture
if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
elseif(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    set(APPLE_ARM TRUE)
elseif(WIN32)
    set(WINDOWS TRUE)
endif()

set(RRTMGP_DATA ${CMAKE_BINARY_DIR}/rrtmgp-data)

# Define the option for enabling/disabling tests
option(ENABLE_TESTS "Enable tests" OFF)

set(RRTMGP_DATA_VERSION "v1.8.2" CACHE STRING "Select what version to clone for rrtmgp-data")
set(RRTMGP_DATA_REPO https://github.com/earth-system-radiation/rrtmgp-data.git)

set(KERNEL_MODE "default" CACHE STRING "Select the kernel mode: default, accel, or extern")
set_property(CACHE KERNEL_MODE PROPERTY STRINGS "default" "accel" "extern")

set(PRECISION "DP" CACHE STRING "Select the floating point precision: double or single")
set_property(CACHE PRECISION PROPERTY STRINGS "DP" "SP")

add_definitions("-DRTE_USE_${PRECISION}")
message(STATUS "Adding compiler definition -DRTE_USE_${PRECISION}")

set(BOOL_TYPE "C" CACHE STRING "Select the bool type: C, Fortran")
set_property(CACHE BOOL_TYPE PROPERTY STRINGS "C" "Fortran")

add_definitions("-DRTE_USE_${BOOL_TYPE}BOOL")
message(STATUS "Adding compiler definition -DRTE_USE_${BOOL_TYPE}BOOL")

set(CMAKE_Fortran_FLAGS "" CACHE STRING "Select Fortran compile flags")

if(NOT WINDOWS)
    set(PREFERRED_FC_COMPILERS "gfortran" "gfortran-10" "gfortran-11" "gfortran-12" "ifort" "ifx" "nvfortran" "f77")
    set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "Select the Fortran compiler: ${PREFERRED_FC_COMPILERS}")
    set_property(CACHE CMAKE_Fortran_COMPILER PROPERTY STRINGS ${PREFERRED_FC_COMPILERS})
else()
    message(STATUS "Conda Environment: $ENV{CONDA_PREFIX}")
    set(CMAKE_Fortran_COMPILER "$ENV{CONDA_PREFIX}/Library/bin/gfortran.exe" CACHE STRING "Select the Fortran compiler")
    set(CMAKE_AR "$ENV{CONDA_PREFIX}/Library/x86_64-w64-mingw32/bin/ar.exe" CACHE STRING "Select archiver")
    set(CMAKE_RANLIB "$ENV{CONDA_PREFIX}/Library/x86_64-w64-mingw32/bin/ranlib.exe" CACHE STRING "Select ranlib executable")

    file(TO_CMAKE_PATH ${CMAKE_AR} CMAKE_AR)
    file(TO_CMAKE_PATH ${CMAKE_RANLIB} CMAKE_RANLIB)
    file(TO_CMAKE_PATH ${CMAKE_Fortran_COMPILER} CMAKE_Fortran_COMPILER)
endif()

# Check if the Fortran compiler is available
find_program(FORTRAN_COMPILER_FOUND ${CMAKE_Fortran_COMPILER})
if(NOT FORTRAN_COMPILER_FOUND)
    message(FATAL_ERROR "No suitable Fortran compiler found. Please install one of the following: ${PREFERRED_FC_COMPILERS}")
endif()

enable_language(Fortran)

include(FetchContent)

message(STATUS "Fetching ${RRTMGP_DATA_REPO} version: ${RRTMGP_DATA_VERSION}")

FetchContent_Declare(
    rrtmgp-data
    GIT_REPOSITORY ${RRTMGP_DATA_REPO}
    GIT_TAG ${RRTMGP_DATA_VERSION}
    SOURCE_DIR ${RRTMGP_DATA}
)
FetchContent_GetProperties(rrtmgp-data)

if(NOT rrtmgp-data_POPULATED)
    FetchContent_MakeAvailable(rrtmgp-data)
endif()

message(STATUS "KERNEL_MODE: ${KERNEL_MODE}")
message(STATUS "PRECISION: ${PRECISION} (DP = 64-bits, SP = 32-bits float)")
message(STATUS "FORTRAN COMPILER: ${CMAKE_Fortran_COMPILER}")

# Set the directory for Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_subdirectory(rte-kernels)     # Build library: rtekernels
add_subdirectory(rte-frontend)    # Build library: rte
add_subdirectory(rrtmgp-kernels)  # Build library: rrtmgpkernels
add_subdirectory(rrtmgp-frontend) # Build library: rrtmgp

if(ENABLE_TESTS)
    find_package(Python REQUIRED COMPONENTS Interpreter)

    # Define failure threshold as an environment variable
    set(FAILURE_THRESHOLD 7.e-4 CACHE STRING "Default failure threshold for tests")

    # Enable testing with CTest
    enable_testing()

    message(STATUS "Tests enabled")

    # Manually locate the directory containing netcdf.mod
    # and the NetCDF Fortran library
    find_path(NETCDF_Fortran_INCLUDE_DIR netcdf.mod HINTS /usr/include $ENV{NFHOME}/include $ENV{CONDA_PREFIX}/include)
    find_library(NETCDF_Fortran_LIBRARY NAMES netcdff HINTS /usr/lib $ENV{NFHOME}/lib $ENV{CONDA_PREFIX}/lib)

    if(NOT NETCDF_Fortran_INCLUDE_DIR OR NOT NETCDF_Fortran_LIBRARY)
        message(SEND_ERROR "NetCDF Fortran not found.")
        if(LINUX AND NOT DEFINED ENV{CONDA_PREFIX})
            message(FATAL_ERROR "You can install it using: `sudo apt-get install libnetcdff-dev`")
        elseif(APPLE_ARM)
            message(FATAL_ERROR "You can install it using: `brew install netcdf`")
        else()
            message(FATAL_ERROR "You can install it using: `conda install -c conda-forge netcdf-fortran`")
        endif()
    else()
        message(STATUS "Found NetCDF module include dir: ${NETCDF_Fortran_INCLUDE_DIR}")
        message(STATUS "Found NetCDF module library dir: ${NETCDF_Fortran_LIBRARY}")
    endif()

    # Function to define a custom test for CTest
    # Arguments:
    #   TEST_NAME - The name of the test in CTest
    #   ARGN      - The libraries that are going to be linked
    function(build_test TEST_NAME)
        set(TEST_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.F90)
        set(TEST_LIBS ${ARGN})

        add_executable(${TEST_NAME} ${TEST_SOURCE})
        target_link_libraries(${TEST_NAME} PRIVATE ${TEST_LIBS})

        message(STATUS "Configuring test ${TEST_NAME}")
    endfunction()

    # Function to define a custom test for CTest
    # Arguments:
    #   TEST_NAME    - The name of the test in CTest
    #   COMMAND_NAME - The command to execute the test
    function(add_custom_test TEST_NAME COMMAND_NAME)
        # Add the test to CTest with a specific command and working directory
        add_test(NAME ${TEST_NAME}
            COMMAND ${COMMAND_NAME} ${ARGN}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )

        message(STATUS "Added test ${TEST_NAME}")
    endfunction()

    set(TEST_ATMOSPHERES ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc)

    # Add the tests subdirectory for building and running tests
    add_subdirectory(examples) # Build: Test & Check
    add_subdirectory(tests)    # Build: Test & Check
endif()
