set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(
  testing-utils STATIC # cmake-format: sort
  mo_simple_netcdf.F90 mo_testing_utils.F90
)

target_include_directories(
  testing-utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)

target_link_libraries(
  testing-utils
  PUBLIC rte-rrtmgp::rrtmgp rte-rrtmgp::rte ${NetCDF_Fortran_LIBRARY}
)

add_library(
  optics-utils STATIC # cmake-format: sort
  mo_optics_utils_rrtmgp.F90
)

target_include_directories(
  optics-utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)

target_link_libraries(
  optics-utils
  PUBLIC rte-rrtmgp::rrtmgp
         rte-rrtmgp::rte
         testing-utils
         ${NetCDF_Fortran_LIBRARY}
)
