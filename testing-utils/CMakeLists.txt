set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

if(TEST_WITH EQUAL "RRTMGP")
  add_library(
    testing-utils STATIC # cmake-format: sort
    mo_optics_utils_rrtmgp.F90 mo_simple_netcdf.F90 mo_testing_utils.F90
  )

  target_include_directories(
    testing-utils
    PUBLIC
      $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
      ${NetCDF_Fortran_INCLUDE_DIR}
  )

elseif(TEST_WITH EQUAL "SSM")
  add_library(
    testing-utils STATIC # cmake-format: sort
    mo_optics_utils_ssm.F90 mo_simple_netcdf.F90 mo_testing_utils.F90
  )

  target_include_directories(
    testing-utils
    PUBLIC
      $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
  )

endif()
