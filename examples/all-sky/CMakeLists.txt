message(STATUS "Configuring lib all_sky_utils")

set(TEST_INPUTS
    ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-lw.nc
    ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-sw.nc
    ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-lw-no-aerosols.nc
    ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-sw-no-aerosols.nc
)

# Copy required test inputs to the build directory
add_custom_target(all_sky_test_inputs ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${TEST_INPUTS} ${CMAKE_CURRENT_BINARY_DIR}/
    COMMENT "Copying required .nc test files to build directory"
)

# Add executable and additional sources
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_load_aerosol_coefficients.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_load_cloud_coefficients.F90
)

add_library(all_sky_utils STATIC ${SOURCES})

target_link_libraries(all_sky_utils PRIVATE
    examples_utils
)

message(STATUS "Configuring lib all_sky_utils")

# Build tests
build_test(rrtmgp_allsky all_sky_utils examples_utils)

add_dependencies(rrtmgp_allsky all_sky_test_inputs)

# Add aerosols tests
add_custom_test(allsky_test_lw rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc)
add_custom_test(allsky_test_sw rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc)

# Add no_aerosols test
add_custom_test(allsky_test_lw_no_aerosols rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc)
add_custom_test(allsky_test_sw_no_aerosols rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw-no-aerosols.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc)

# Add aerosols and no_aerosols check test
add_custom_test(check_allsky_lw_sw
    ${CMAKE_COMMAND} -E env ${Python_EXECUTABLE}
        ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
        --ref_dir ${RRTMGP_DATA}/examples/all-sky/reference
        --tst_dir ${CMAKE_CURRENT_BINARY_DIR}
        --variables lw_flux_up lw_flux_dn sw_flux_up sw_flux_dn sw_flux_dir
        --file_names rrtmgp-allsky-lw.nc rrtmgp-allsky-sw.nc
        --failure_threshold ${FAILURE_THRESHOLD}
)
add_custom_test(check_allsky_no_aerosols_lw_sw
    ${CMAKE_COMMAND} -E env ${Python_EXECUTABLE}
        ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
        --ref_dir ${RRTMGP_DATA}/examples/all-sky/reference
        --tst_dir ${CMAKE_CURRENT_BINARY_DIR}
        --variables lw_flux_up lw_flux_dn sw_flux_up sw_flux_dn sw_flux_dir
        --file_names rrtmgp-allsky-lw-no-aerosols.nc rrtmgp-allsky-sw-no-aerosols.nc
        --failure_threshold ${FAILURE_THRESHOLD}
)
