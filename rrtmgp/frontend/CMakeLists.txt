set(gas_optics_source_dir ${PROJECT_SOURCE_DIR}/gas-optics)

add_library(
  rrtmgp # cmake-format: sort
  $<TARGET_OBJECTS:rrtmgpkernels>
  ${gas_optics_source_dir}/mo_gas_concentrations.F90
  ${gas_optics_source_dir}/mo_gas_optics.F90
  ${gas_optics_source_dir}/mo_gas_optics_constants.F90
  ${gas_optics_source_dir}/mo_gas_optics_util_string.F90
  mo_aerosol_optics_rrtmgp_merra.F90
  mo_cloud_optics_rrtmgp.F90
  mo_gas_optics_rrtmgp.F90
)

add_library(rte-rrtmgp::rrtmgp ALIAS rrtmgp)

set_target_properties(
  rrtmgp
  PROPERTIES EXPORT_NAME rte-rrtmgp::rrtmgp
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
             INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}"
)

target_include_directories(
  rrtmgp
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
)

target_link_libraries(rrtmgp PRIVATE rte)

install(TARGETS rrtmgp EXPORT rte-rrtmgp-targets)
