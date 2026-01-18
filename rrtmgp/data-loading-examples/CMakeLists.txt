set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

add_library(
  rrtmgp-load-data STATIC # cmake-format: sort
  mo_optics_utils_rrtmgp.F90
)
target_link_libraries(
  rrtmgp-load-data PUBLIC shared-utils rte-rrtmgp::rrtmgp rte-rrtmgp::rte
)
target_include_directories(
  rrtmgp-load-data
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)
