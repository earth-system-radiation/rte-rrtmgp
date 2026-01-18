set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(
  shared-utils STATIC # cmake-format: sort
  mo_simple_netcdf.F90 mo_testing_utils.F90
)
target_link_libraries(
  shared-utils
  PUBLIC rte-rrtmgp::rrtmgp rte-rrtmgp::rte ${NetCDF_Fortran_LIBRARY}
)
target_include_directories(
  shared-utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)
