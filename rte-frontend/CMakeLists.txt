set(RTE_KIND_SRC
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_kind.F90
)

SET(SOURCES
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_fluxes.F90
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_optical_props.F90
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_config.F90
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_lw.F90
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_sw.F90
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_util_array_validation.F90
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_source_functions.F90
)

add_library(rtef_kind STATIC ${RTE_KIND_SRC})

message("-- Building rtef_kind.a")

add_library(rtef STATIC ${SOURCES})
target_link_libraries(rtef
    PRIVATE
        rtef_kind
)

message("-- Building rtef.a")
