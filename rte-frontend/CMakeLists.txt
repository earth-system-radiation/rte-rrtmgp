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

add_library(rtef_objects OBJECT ${SOURCES})
add_library(rtekind_object OBJECT ${RTE_KIND_SRC})

add_library(rtef STATIC $<TARGET_OBJECTS:rtekind_object> $<TARGET_OBJECTS:rtef_objects>)

message("-- Building rtef.a")
