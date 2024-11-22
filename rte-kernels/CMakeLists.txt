SET(SOURCES
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_rte_solver_kernels.F90
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_fluxes_broadband_kernels.F90
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_rte_util_array.F90
    ${CMAKE_SOURCE_DIR}/rte-kernels/mo_optical_props_kernels.F90
)

add_library(rtekernels STATIC ${SOURCES})

message("-- Building rtekernels.a")
