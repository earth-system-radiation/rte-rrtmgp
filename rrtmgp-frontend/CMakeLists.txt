SET(SOURCES
    ${CMAKE_SOURCE_DIR}/gas-optics/mo_gas_concentrations.F90
    ${CMAKE_SOURCE_DIR}/gas-optics/mo_gas_optics_constants.F90
    ${CMAKE_SOURCE_DIR}/gas-optics/mo_gas_optics_util_string.F90
    ${CMAKE_SOURCE_DIR}/gas-optics/mo_gas_optics.F90

    ${CMAKE_SOURCE_DIR}/rrtmgp-frontend/mo_aerosol_optics_rrtmgp_merra.F90
    ${CMAKE_SOURCE_DIR}/rrtmgp-frontend/mo_cloud_optics_rrtmgp.F90
    ${CMAKE_SOURCE_DIR}/rrtmgp-frontend/mo_gas_optics_rrtmgp.F90
)

add_library(rrtmgpf_objects OBJECT ${SOURCES})
add_library(rrtmgpf STATIC $<TARGET_OBJECTS:rrtmgpf_objects>)

message("-- Building rrtmgpf.a")
