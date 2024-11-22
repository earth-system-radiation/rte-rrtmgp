SET(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_rte_config.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_rte_util_array_validation.F90

    ${CMAKE_CURRENT_SOURCE_DIR}/mo_rte_lw.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_rte_sw.F90

    # FLUX REDUCTION
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_fluxes.F90
    # OPTICAL PROPERTIES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_optical_props.F90
    # SOURCE FUNCTIONS
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_source_functions.F90
)

add_library(rte STATIC ${SOURCES})
target_link_libraries(rte PUBLIC
    rtef_kind
    rtekernels
)

message(STATUS "Building lib rte.a")
