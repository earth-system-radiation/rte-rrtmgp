set(gas_optics_source_dir ./gas-optics-template)

add_library(
  rte # cmake-format: sort
  $<TARGET_OBJECTS:rtekernels>
  ${gas_optics_source_dir}/mo_gas_concentrations.F90
  ${gas_optics_source_dir}/mo_gas_optics.F90
  ${gas_optics_source_dir}/mo_gas_optics_constants.F90
  ${gas_optics_source_dir}/mo_gas_optics_util_string.F90
  mo_fluxes.F90
  mo_optical_props.F90
  mo_rte_config.F90
  mo_rte_lw.F90
  mo_rte_sw.F90
  mo_rte_util_array_validation.F90
  mo_source_functions.F90
)

add_library(rte-rrtmgp::rte ALIAS rte)

set_target_properties(
  rte
  PROPERTIES EXPORT_NAME rte-rrtmgp::rte
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

target_include_directories(
  rte
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)

install(TARGETS rte EXPORT rte-rrtmgp-targets)
