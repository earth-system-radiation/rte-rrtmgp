set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

if(TEST_WITH EQUAL RRTMGP)
  set(optics-interface mo_optics_utils_rrtmgp.F90)
elseif(TEST_WITH EQUAL SSM)
  set(optics-interface mo_optics_utils_ssm.F90)
endif()

# Need to figure out how to pick the right optics interface from TEST_WITH
set(optics-interface mo_optics_utils_rrtmgp.F90)

add_library(
  testing-utils STATIC # cmake-format: sort
  ${optics-interface} mo_simple_netcdf.F90 mo_testing_utils.F90
)
target_link_libraries(
  testing-utils
  PUBLIC rte-rrtmgp::rrtmgp rte-rrtmgp::rte ${NetCDF_Fortran_LIBRARY}
)
target_include_directories(
  testing-utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)
