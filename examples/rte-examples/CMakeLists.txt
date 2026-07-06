#
# Example-specific utilities. Need to specify include directories because
# mo_rte_examples_io.F90 uses netcdf directly
#
set(executable rte_examples)

add_library(${executable}_io STATIC mo_${executable}_io.F90)
target_include_directories(
  ${executable}_io
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)
target_link_libraries(rte_examples_io PRIVATE shared-utils)

# ##############################################################################
#
# RTE tests using homogenized examples
#
# Build the two executables
#
add_executable(${executable} ${executable}.F90)
target_include_directories(
  ${executable}
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)
target_link_libraries(
  ${executable}
  PRIVATE ${executable}_io
          rrtmgp-load-data
          shared-utils
          rte-rrtmgp::rte
          rte-rrtmgp::rrtmgp
          rte-rrtmgp::ssm
)
