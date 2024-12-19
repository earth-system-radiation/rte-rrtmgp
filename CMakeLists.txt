cmake_minimum_required(VERSION 3.18)

project(rte-rrtmgp LANGUAGES Fortran)

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
