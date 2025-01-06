cmake_minimum_required(VERSION 3.18)

project(
  rte-rrtmgp
  VERSION 1.9.0
  LANGUAGES Fortran
)

option(BUILD_TESTING "Build tests" OFF)

option(RTE_ENABLE_SP "Enable single-precision floating-point model" OFF)

set(PREFERRED_KERNEL_MODES "default" "accel" "extern")
set(KERNEL_MODE
    "default"
    CACHE STRING "Select the kernel mode: ${PREFERRED_KERNEL_MODES}"
)
set_property(CACHE KERNEL_MODE PROPERTY STRINGS ${PREFERRED_KERNEL_MODES})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_compile_options(
  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>
)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

add_subdirectory(rte-kernels)
add_subdirectory(rte-frontend)
add_subdirectory(rrtmgp-kernels)
add_subdirectory(rrtmgp-frontend)

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)

  if(NOT _RTE_RRTMGP_HAVE_PY_PACKAGES)
    include(CheckPython3Package)
    check_python3_package(numpy)
    check_python3_package(
      "netCDF4 or h5netcdf+scipy"
      CODE "try:
    import netCDF4
except:
    import h5netcdf
    import scipy
"
    )
    check_python3_package(
      "xarray>=0.12.2"
      CODE "import xarray
exit(tuple(map(int, xarray.__version__.split('.'))) < (0, 12, 2))"
    )
    check_python3_package(dask.array)
    set(_RTE_RRTMGP_HAVE_PY_PACKAGES
        TRUE
        CACHE INTERNAL
              "RTE-RRTMGP found all Python packages required for testing"
    )
  endif()

  find_package(NetCDF_Fortran REQUIRED)

  if(RRTMGP_DATA)
    add_test(NAME fetch_rrtmgp_data COMMAND ${CMAKE_COMMAND} -E true)
  else()
    set(RRTMGP_DATA "${PROJECT_BINARY_DIR}/rrtmgp-data")

    include(ExternalProject)
    ExternalProject_Add(
      rrtmgp-data
      GIT_REPOSITORY https://github.com/earth-system-radiation/rrtmgp-data.git
      GIT_TAG "v1.8.2"
      GIT_SHALLOW True
      EXCLUDE_FROM_ALL True
      PREFIX rrtmgp-data-cmake
      SOURCE_DIR ${RRTMGP_DATA}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
    )
    add_test(
      NAME fetch_rrtmgp_data
      COMMAND
        ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --config
        "$<CONFIG>" --target rrtmgp-data
    )
  endif()

  set_tests_properties(
    fetch_rrtmgp_data PROPERTIES FIXTURES_SETUP fetch_rrtmgp_data
  )

  add_subdirectory(examples)
  add_subdirectory(tests)
else()
  # Allow for 'make test' even if the tests are disabled:
  enable_testing()
endif()
