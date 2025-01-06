if(RTE_ENABLE_SP)
  set(default_FAILURE_THRESHOLD "3.5e-1")
else()
  set(default_FAILURE_THRESHOLD "7.e-4")
endif()

set(FAILURE_THRESHOLD
    "${default_FAILURE_THRESHOLD}"
    CACHE STRING "Default failure threshold for tests"
)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(
  examples_utils STATIC # cmake-format: sort
  mo_garand_atmos_io.F90 mo_load_coefficients.F90 mo_simple_netcdf.F90
)

target_include_directories(
  examples_utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)

target_link_libraries(
  examples_utils
  PUBLIC rte-rrtmgp::rrtmgp rte-rrtmgp::rte ${NetCDF_Fortran_LIBRARY}
)

add_subdirectory(all-sky)
add_subdirectory(rfmip-clear-sky)
