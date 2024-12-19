cmake_minimum_required(VERSION 3.18)

project(
  rte-rrtmgp
  VERSION 1.8
  LANGUAGES Fortran
)

option(BUILD_TESTING "Build tests" OFF)

option(RTE_ENABLE_SP "Enable single-precision floating-point model" OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
)

set(PREFERRED_KERNEL_MODES "default" "accel" "extern")
set(KERNEL_MODE
    "default"
    CACHE STRING "Select the kernel mode: ${PREFERRED_KERNEL_MODES}"
)
set_property(CACHE KERNEL_MODE PROPERTY STRINGS ${PREFERRED_KERNEL_MODES})

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

add_subdirectory(rte-kernels)
add_subdirectory(rte-frontend)
add_subdirectory(rrtmgp-kernels)
add_subdirectory(rrtmgp-frontend)

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  find_package(NetCDF_Fortran REQUIRED)

  set(RRTMGP_DATA ${PROJECT_BINARY_DIR}/rrtmgp-data)
  set(RRTMGP_DATA_VERSION "v1.8.2")
  include(FetchContent)
  message(CHECK_START "Fetching RRTMGP data ${RRTMGP_DATA_VERSION}")
  FetchContent_Declare(
    rrtmgp-data
    GIT_REPOSITORY https://github.com/earth-system-radiation/rrtmgp-data.git
    GIT_TAG ${RRTMGP_DATA_VERSION}
    SOURCE_DIR ${RRTMGP_DATA}
  )
  FetchContent_MakeAvailable(rrtmgp-data)
  message(CHECK_PASS "done")
  add_subdirectory(examples)
  add_subdirectory(tests)
else()
  # Allow for 'make test' even if the tests are disabled:
  enable_testing()
endif()
