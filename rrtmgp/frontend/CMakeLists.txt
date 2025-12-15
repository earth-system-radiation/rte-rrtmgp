add_library(
  rrtmgp # cmake-format: sort
  $<TARGET_OBJECTS:rrtmgpkernels>
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
