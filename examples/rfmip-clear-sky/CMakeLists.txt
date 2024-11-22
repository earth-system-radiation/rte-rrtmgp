set(TEST_INPUTS
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
)

# Copy required test inputs to the build directory
add_custom_target(rfmip_clear_sky_test_inputs ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${TEST_INPUTS} ${CMAKE_CURRENT_BINARY_DIR}/
    COMMENT "Copying required .nc test files to build directory"
)

# Add executable and additional sources
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_rfmip_io.F90
)

add_library(rfmip_clear_utils STATIC ${SOURCES})

target_link_libraries(rfmip_clear_utils PRIVATE
    examples_utils
)

message(STATUS "Building lib rfmip_clear_utils.a")

# Path to the test atmospheres file

# Build tests
build_test(rrtmgp_rfmip_lw rfmip_clear_utils examples_utils)
build_test(rrtmgp_rfmip_sw rfmip_clear_utils examples_utils)

add_dependencies(rrtmgp_rfmip_lw rfmip_clear_sky_test_inputs)
add_dependencies(rrtmgp_rfmip_sw rfmip_clear_sky_test_inputs)

# Add rfmip tests
add_custom_test(rfmip_test_lw rrtmgp_rfmip_lw 8 multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc)
add_custom_test(rfmip_test_sw rrtmgp_rfmip_sw 8 multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc)

# Add rfmip check test
add_custom_test(check_rfmip_lw_sw
    ${CMAKE_COMMAND} -E env ${Python_EXECUTABLE}
        ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
        --ref_dir ${RRTMGP_DATA}/examples/rfmip-clear-sky/reference
        --tst_dir ${CMAKE_CURRENT_BINARY_DIR}
        --variables rld rlu rsd rsu
        --file_names rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
        --failure_threshold ${FAILURE_THRESHOLD}
)
