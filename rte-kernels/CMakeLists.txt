SET(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (RTE_KERNELS STREQUAL "ACCEL")
    SET(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/accel)
elseif(RTE_KERNELS STREQUAL "EXTERN")
    SET(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/api)
endif()

if (RTE_KERNELS STREQUAL "ACCEL")
    SET(FLUX_REDUCTION  ${CMAKE_CURRENT_SOURCE_DIR}/mo_fluxes_broadband_kernels.F90)
    SET(ARRAY_UTILITIES ${CMAKE_CURRENT_SOURCE_DIR}/mo_rte_util_array.F90)
else()
    SET(FLUX_REDUCTION  ${LOCAL_SOURCE_DIR}/mo_fluxes_broadband_kernels.F90)
    SET(ARRAY_UTILITIES ${LOCAL_SOURCE_DIR}/mo_rte_util_array.F90)
endif()

SET(OPTICAL_PROPERTIES ${LOCAL_SOURCE_DIR}/mo_optical_props_kernels.F90)
SET(RADIATIVE_TRANSFER ${LOCAL_SOURCE_DIR}/mo_rte_solver_kernels.F90)

SET(SOURCES
    ${FLUX_REDUCTION}
    ${ARRAY_UTILITIES}
    # OPTICAL PROPERTIES
    ${LOCAL_SOURCE_DIR}/mo_optical_props_kernels.F90
    # RADIATIVE TRANSFER
    ${LOCAL_SOURCE_DIR}/mo_rte_solver_kernels.F90
)

add_library(rtekernels STATIC ${SOURCES})
target_link_libraries(rtekernels PRIVATE
    rtef_kind
)

message(STATUS "Building lib rtekernels.a")
