cmake_minimum_required(VERSION 3.4...3.18)
project(rte-rrtmgp LANGUAGES)

# Define flags to identify the operating system
set(LINUX FALSE)
set(APPLE_ARM FALSE)
set(WINDOWS FALSE)

# Determine the operating system and architecture
if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
elseif(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
  set(APPLE_ARM TRUE)
elseif(WIN32)
  set(WINDOWS TRUE)
endif()

# Check if CMAKE_BUILD_TYPE is defined; if not, set a default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING
            "Build type: Debug, Release, RelWithDebInfo, or MinSizeRel"
            FORCE
  )
endif()

# Provide the list of valid build types
set_property(
  CACHE CMAKE_BUILD_TYPE
  PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel"
)

# Inform the user about the selected build type
message(STATUS "Using build type: ${CMAKE_BUILD_TYPE}")

set(RRTMGP_DATA ${CMAKE_BINARY_DIR}/rrtmgp-data)

# Define the option for enabling/disabling tests
option(ENABLE_TESTS "Enable tests" OFF)

set(RRTMGP_DATA_VERSION
    "v1.8.2"
    CACHE STRING "Select what version to clone for rrtmgp-data"
)
set(RRTMGP_DATA_REPO https://github.com/earth-system-radiation/rrtmgp-data.git)

set(PREFERRED_KERNEL_MODES "default" "accel" "extern")

set(KERNEL_MODE
    "default"
    CACHE STRING "Select the kernel mode: ${PREFERRED_KERNEL_MODES}"
)
set_property(CACHE KERNEL_MODE PROPERTY STRINGS ${PREFERRED_KERNEL_MODES})

set(PRECISION
    "DP"
    CACHE STRING "Select the floating point precision: double or single"
)
set_property(CACHE PRECISION PROPERTY STRINGS "DP" "SP")

add_definitions("-DRTE_USE_${PRECISION}")
message(STATUS "Adding compiler definition -DRTE_USE_${PRECISION}")

option(USE_C_BOOL "Use C bool" OFF)

if(USE_C_BOOL)
  add_definitions("-DRTE_USE_CBOOL")
  message(STATUS "Adding compiler definition -DRTE_USE_CBOOL")
endif()

set(CMAKE_Fortran_FLAGS
    ""
    CACHE STRING "Select Fortran compile flags"
)

set(PREFERRED_FC_COMPILERS
    "gfortran"
    "gfortran-10"
    "gfortran-11"
    "gfortran-12"
    "gfortran-13"
    "ifort"
    "ifx"
    "nvfortran"
    "f77"
)

if(NOT WINDOWS)
  set(CMAKE_Fortran_COMPILER
      "gfortran"
      CACHE STRING "Select the Fortran compiler: ${PREFERRED_FC_COMPILERS}"
  )
  set_property(
    CACHE CMAKE_Fortran_COMPILER PROPERTY STRINGS ${PREFERRED_FC_COMPILERS}
  )
else()
  message(STATUS "Conda Environment: $ENV{CONDA_PREFIX}")
  set(CMAKE_Fortran_COMPILER
      "$ENV{CONDA_PREFIX}/Library/bin/gfortran.exe"
      CACHE STRING "Select the Fortran compiler"
  )
  set(CMAKE_AR
      "$ENV{CONDA_PREFIX}/Library/x86_64-w64-mingw32/bin/ar.exe"
      CACHE STRING "Select archiver"
  )
  set(CMAKE_RANLIB
      "$ENV{CONDA_PREFIX}/Library/x86_64-w64-mingw32/bin/ranlib.exe"
      CACHE STRING "Select ranlib executable"
  )

  file(TO_CMAKE_PATH ${CMAKE_AR} CMAKE_AR)
  file(TO_CMAKE_PATH ${CMAKE_RANLIB} CMAKE_RANLIB)
  file(TO_CMAKE_PATH ${CMAKE_Fortran_COMPILER} CMAKE_Fortran_COMPILER)
endif()

# Check if the Fortran compiler is available
find_program(FORTRAN_COMPILER_FOUND ${CMAKE_Fortran_COMPILER})
if(NOT FORTRAN_COMPILER_FOUND)
  message(
    FATAL_ERROR
      "No suitable Fortran compiler found. Please install one of the following: ${PREFERRED_FC_COMPILERS}"
  )
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

message(STATUS "KERNEL_MODE: ${KERNEL_MODE} (${PREFERRED_KERNEL_MODES})")
message(STATUS "PRECISION: ${PRECISION} (DP = 64-bits, SP = 32-bits float)")
message(STATUS "FORTRAN COMPILER: ${CMAKE_Fortran_COMPILER}")

# Set the directory for Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

include(common.cmake)

add_subdirectory(rte-kernels) # Build library: rtekernels
add_subdirectory(rte-frontend) # Build library: rte
add_subdirectory(rrtmgp-kernels) # Build library: rrtmgpkernels
add_subdirectory(rrtmgp-frontend) # Build library: rrtmgp

if(ENABLE_TESTS)
  include(tests.cmake)

  # Add the tests subdirectory for building and running tests
  add_subdirectory(examples) # Build: Test & Check
  add_subdirectory(tests) # Build: Test & Check
endif()
