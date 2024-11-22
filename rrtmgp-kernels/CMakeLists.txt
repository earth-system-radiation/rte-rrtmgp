SET(SOURCES
    # "../rte-frontend/mo_rte_kind.F90"
    ${CMAKE_SOURCE_DIR}/rrtmgp-kernels/mo_gas_optics_rrtmgp_kernels.F90
)

add_library(rrtmgpkernels_objects OBJECT ${SOURCES})
add_library(rrtmgpkernels STATIC $<TARGET_OBJECTS:rrtmgpkernels_objects>)

message("-- Building rrtmgpkernels.a")
