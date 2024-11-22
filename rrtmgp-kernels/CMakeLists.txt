SET(LOCAL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (RTE_KERNELS STREQUAL "ACCEL")
    SET(LOCAL_SOURCE_DIR ${LOCAL_SOURCE_DIR}/accel)
elseif(RTE_KERNELS STREQUAL "EXTERN")
    SET(LOCAL_SOURCE_DIR ${LOCAL_SOURCE_DIR}/api)
endif()

SET(SOURCES
    ${LOCAL_SOURCE_DIR}/mo_gas_optics_rrtmgp_kernels.F90
)

add_library(rrtmgpkernels STATIC ${SOURCES})
target_include_directories(rrtmgpkernels PUBLIC
    ${LOCAL_SOURCE_DIR}
)
target_link_libraries(rrtmgpkernels PRIVATE
    rtekernels
)

message(STATUS "Building lib rrtmgpkernels.a")
