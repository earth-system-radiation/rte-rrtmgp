add_library(
  rte STATIC # cmake-format: sort
  mo_fluxes.F90
  mo_optical_props.F90
  mo_rte_config.F90
  mo_rte_lw.F90
  mo_rte_sw.F90
  mo_rte_util_array_validation.F90
  mo_source_functions.F90
)

target_include_directories(
  rte
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)

target_link_libraries(rte PUBLIC rtekernels)
