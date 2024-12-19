set(gas_optics_source_dir ${PROJECT_SOURCE_DIR}/gas-optics)

add_library(
  rrtmgp STATIC # cmake-format: sort
  ${gas_optics_source_dir}/mo_gas_concentrations.F90
  ${gas_optics_source_dir}/mo_gas_optics.F90
  ${gas_optics_source_dir}/mo_gas_optics_constants.F90
  ${gas_optics_source_dir}/mo_gas_optics_util_string.F90
  mo_aerosol_optics_rrtmgp_merra.F90
  mo_cloud_optics_rrtmgp.F90
  mo_gas_optics_rrtmgp.F90
)

target_link_libraries(rrtmgp PUBLIC rrtmgpkernels rte)
