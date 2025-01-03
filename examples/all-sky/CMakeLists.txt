add_library(
  all_sky_utils STATIC # cmake-format: sort
  mo_load_aerosol_coefficients.F90 mo_load_cloud_coefficients.F90
)

target_link_libraries(all_sky_utils PUBLIC examples_utils)

add_executable(rrtmgp_allsky rrtmgp_allsky.F90)
target_link_libraries(rrtmgp_allsky PRIVATE all_sky_utils)

add_test(
  NAME run_allsky_lw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc
    ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc
)
set_tests_properties(
  run_allsky_lw
  PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data FIXTURES_SETUP run_allsky
)

add_test(
  NAME run_allsky_sw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc
    ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc
)
set_tests_properties(
  run_allsky_sw
  PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data FIXTURES_SETUP run_allsky
)

add_test(
  NAME check_allsky_lw_sw
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
    --ref_dir ${RRTMGP_DATA}/examples/all-sky/reference --tst_dir
    ${CMAKE_CURRENT_BINARY_DIR} --variables lw_flux_up lw_flux_dn sw_flux_up
    sw_flux_dn sw_flux_dir --file_names rrtmgp-allsky-lw.nc rrtmgp-allsky-sw.nc
    --failure_threshold ${FAILURE_THRESHOLD}
)
set_tests_properties(
  check_allsky_lw_sw
  PROPERTIES FIXTURES_REQUIRED "fetch_rrtmgp_data;run_allsky"
)

add_test(
  NAME run_allsky_no_aerosols_lw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc
)
set_tests_properties(
  run_allsky_no_aerosols_lw
  PROPERTIES FIXTURES_REQUIRED
             fetch_rrtmgp_data
             FIXTURES_SETUP
             run_allsky_no_aerosols
)

add_test(
  NAME run_allsky_no_aerosols_sw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw-no-aerosols.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc
)
set_tests_properties(
  run_allsky_no_aerosols_sw
  PROPERTIES FIXTURES_REQUIRED
             fetch_rrtmgp_data
             FIXTURES_SETUP
             run_allsky_no_aerosols
)

add_test(
  NAME check_allsky_no_aerosols_lw_sw
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
    --ref_dir ${RRTMGP_DATA}/examples/all-sky/reference --tst_dir
    ${CMAKE_CURRENT_BINARY_DIR} --variables lw_flux_up lw_flux_dn sw_flux_up
    sw_flux_dn sw_flux_dir --file_names rrtmgp-allsky-lw-no-aerosols.nc
    rrtmgp-allsky-sw-no-aerosols.nc --failure_threshold ${FAILURE_THRESHOLD}
)
set_tests_properties(
  check_allsky_no_aerosols_lw_sw
  PROPERTIES FIXTURES_REQUIRED "fetch_rrtmgp_data;run_allsky_no_aerosols"
)
