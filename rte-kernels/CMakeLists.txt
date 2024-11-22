SET(SOURCES
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_rte_solver_kernels.F90
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_fluxes_broadband_kernels.F90
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_rte_util_array.F90
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_optical_props_kernels.F90
)

add_library(rtekernels_objects OBJECT ${SOURCES})
add_library(rtekernels STATIC $<TARGET_OBJECTS:rtekind_object> $<TARGET_OBJECTS:rtekernels_objects>)

message("-- Building rtekernels.a")
