SET(GAS_OPTICS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/gas-optics)

SET(SOURCES
    # PHYSICAL CONSTANTS
    ${GAS_OPTICS_SOURCE_DIR}/mo_gas_optics_constants.F90
    # UTILITY
    ${GAS_OPTICS_SOURCE_DIR}/mo_gas_optics_util_string.F90
    # GAS CONCENTRATIONS
    ${GAS_OPTICS_SOURCE_DIR}/mo_gas_concentrations.F90
    # GAS OPTICS
    ${GAS_OPTICS_SOURCE_DIR}/mo_gas_optics.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_gas_optics_rrtmgp.F90
)

add_library(rrtmgp STATIC ${SOURCES})
target_link_libraries(rrtmgp PRIVATE
    rrtmgpkernels
    rte
)

message(STATUS "Building lib rrtmgp.a")
