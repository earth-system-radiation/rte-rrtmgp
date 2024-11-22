set(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (RTE_KERNELS STREQUAL "accel")
    set(LOCAL_SOURCE_DIR ${LOCAL_SOURCE_DIR}/accel)
elseif(RTE_KERNELS STREQUAL "extern")
    set(LOCAL_SOURCE_DIR ${LOCAL_SOURCE_DIR}/api)
endif()

set(SOURCES
    ${LOCAL_SOURCE_DIR}/mo_gas_optics_rrtmgp_kernels.F90
)

add_library(rrtmgpkernels STATIC ${SOURCES})
target_link_libraries(rrtmgpkernels PRIVATE
    rtekernels
)

message(STATUS "Building lib rrtmgpkernels.a")
