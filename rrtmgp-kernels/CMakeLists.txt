set(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (KERNEL_MODE STREQUAL "accel")
    set(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/accel)
elseif(KERNEL_MODE STREQUAL "extern")
    set(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/api)
endif()

if (KERNEL_MODE STREQUAL "accel")
    set(CLOUD_OPTICS ${CMAKE_CURRENT_SOURCE_DIR}/mo_cloud_optics_rrtmgp_kernels.F90)
else()
    set(CLOUD_OPTICS ${LOCAL_SOURCE_DIR}/mo_cloud_optics_rrtmgp_kernels.F90)
endif()

set(GAS_OPTICS ${LOCAL_SOURCE_DIR}/mo_gas_optics_rrtmgp_kernels.F90)

set(SOURCES
    ${CLOUD_OPTICS}
    ${GAS_OPTICS}
)

add_library(rrtmgpkernels STATIC ${SOURCES})
target_link_libraries(rrtmgpkernels PRIVATE
    rtekernels
)

message(STATUS "Building lib rrtmgpkernels.a")
