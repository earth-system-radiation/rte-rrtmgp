cmake_minimum_required(VERSION 3.18)

project(
  rte-rrtmgp
  VERSION 1.8
  LANGUAGES Fortran
)

set(RRTMGP_DATA ${PROJECT_BINARY_DIR}/rrtmgp-data)

option(BUILD_TESTING "Build tests" OFF)

option(RTE_ENABLE_SP "Enable single-precision floating-point model" OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
)

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

if(RTE_ENABLE_SP)
  add_definitions("-DRTE_USE_SP")
endif()

option(USE_C_BOOL "Use C bool" OFF)

if(USE_C_BOOL)
  add_definitions("-DRTE_USE_CBOOL")
  message(STATUS "Adding compiler definition -DRTE_USE_CBOOL")
endif()

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

# Set the directory for Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

add_subdirectory(rte-kernels)
add_subdirectory(rte-frontend)
add_subdirectory(rrtmgp-kernels)
add_subdirectory(rrtmgp-frontend)

include(CTest)
if(BUILD_TESTING)
  find_package(Python REQUIRED COMPONENTS Interpreter)
  find_package(NetCDF_Fortran REQUIRED)
  set(FAILURE_THRESHOLD
      7.e-4
      CACHE STRING "Default failure threshold for tests"
  )
  set(TEST_ATMOSPHERES
      ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  )
  add_subdirectory(examples)
  add_subdirectory(tests)
else()
  # Allow for 'make test' even if the tests are disabled:
  enable_testing()
endif()
