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

set(RTE_KERNEL_MODE "default" CACHE STRING "Select the kernel mode: default, accel, or extern")
set_property(CACHE RTE_KERNEL_MODE PROPERTY STRINGS "default" "accel" "extern")

if("${FP_MODEL}" STREQUAL "SP")
    set(PRECISION "SP" CACHE STRING "Select the floating point precision: double or single" FORCE)
elseif("${FP_MODEL}" STREQUAL "DP")
    set(PRECISION "DP" CACHE STRING "Select the floating point precision: double or single" FORCE)
else()
    set(PRECISION "DP" CACHE STRING "Select the floating point precision: double or single")
endif()
set_property(CACHE PRECISION PROPERTY STRINGS "DP" "SP")

add_definitions("-DRTE_USE_${PRECISION}")
message(STATUS "Adding compiler definition -DRTE_USE_${PRECISION}")

set(CMAKE_Fortran_FLAGS "" CACHE STRING "Select Fortran compile flags")

if(NOT WINDOWS)
    set(PREFERRED_FC_COMPILERS "gfortran" "ifort" "gfortran-10" "gfortran-11" "gfortran-12" "f77")
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

FetchContent_Declare(
    rrtmgp-data
    GIT_REPOSITORY https://github.com/earth-system-radiation/rrtmgp-data.git
    GIT_TAG origin/main
    SOURCE_DIR ${RRTMGP_DATA}
)
FetchContent_GetProperties(rrtmgp-data)

if(NOT rrtmgp-data_POPULATED)
    FetchContent_MakeAvailable(rrtmgp-data)
endif()

message(STATUS "RTE_KERNEL_MODE: ${RTE_KERNEL_MODE}")
message(STATUS "PRECISION: ${PRECISION} (DP = 64-bits, SP = 32-bits float)")
message(STATUS "FORTRAN COMPILER: ${CMAKE_Fortran_COMPILER}")

# Set the directory for Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(RTE_KIND_SRC ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_kind.F90)
add_library(rtef_kind STATIC ${RTE_KIND_SRC})

message(STATUS "Building lib rtef_kind.a")

add_subdirectory(rte-kernels)     # Build library: rtekernels.a
add_subdirectory(rte-frontend)    # Build library: rte.a
add_subdirectory(rrtmgp-kernels)  # Build library: rrtmgpkernels.a
add_subdirectory(rrtmgp-frontend) # Build library: rrtmgp.a

if(ENABLE_TESTS)
    # Enable testing with CTest
    enable_testing()

    message(STATUS "Tests enabled")

    # Manually locate the directory containing netcdf.mod
    # and the NetCDF Fortran library
    if (NOT WIN32)
        find_path(NETCDF_Fortran_INCLUDE_DIR netcdf.mod HINTS /usr/include)
        find_library(NETCDF_Fortran_LIBRARY NAMES netcdff HINTS /usr/lib)
    else()
        find_path(NETCDF_Fortran_INCLUDE_DIR netcdf.mod HINTS $ENV{CONDA_PREFIX}/include)
        find_library(NETCDF_Fortran_LIBRARY NAMES netcdff HINTS $ENV{CONDA_PREFIX}/lib)
    endif()

    if(NOT NETCDF_Fortran_INCLUDE_DIR OR NOT NETCDF_Fortran_LIBRARY)
        message(SEND_ERROR "NetCDF Fortran not found.")
        if(LINUX)
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

        message(STATUS "Building test ${TEST_NAME}.so")
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

    # Add the tests subdirectory for building and running tests
    add_subdirectory(examples) # Build: Test & Check
    add_subdirectory(tests)    # Build: Test & Check
endif()
