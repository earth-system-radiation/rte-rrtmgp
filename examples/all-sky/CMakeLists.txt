# Add executable and additional sources
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_load_aerosol_coefficients.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_load_cloud_coefficients.F90
)

add_library(all_sky_utils STATIC ${SOURCES})

target_link_libraries(all_sky_utils PRIVATE
    examples_utils
)

message(STATUS "Building lib all_sky_utils.a")

set(TEST_ALLSKY_LW             ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-lw.nc)
set(TEST_ALLSKY_SW             ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-sw.nc)
set(TEST_ALLSKY_LW_NO_AEROSOLS ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-lw-no-aerosols.nc)
set(TEST_ALLSKY_SW_NO_AEROSOLS ${RRTMGP_DATA}/examples/all-sky/reference/rrtmgp-allsky-sw-no-aerosols.nc)

# Build tests
build_test(rrtmgp_allsky all_sky_utils examples_utils)

# Add aerosols tests
add_custom_test(rrtmgp_allsky_test_lw rrtmgp_allsky 24 72 1 ${TEST_ALLSKY_LW} ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc)
add_custom_test(rrtmgp_allsky_test_sw rrtmgp_allsky 24 72 1 ${TEST_ALLSKY_SW} ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc)

# Add no_aerosols test
add_custom_test(rrtmgp_allsky_test_lw_no_aerosols rrtmgp_allsky 24 72 1 ${TEST_ALLSKY_LW_NO_AEROSOLS} ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc)
add_custom_test(rrtmgp_allsky_test_sw_no_aerosols rrtmgp_allsky 24 72 1 ${TEST_ALLSKY_SW_NO_AEROSOLS} ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc)
