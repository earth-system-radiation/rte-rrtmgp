set(rrtmgp_load_source_dir ${PROJECT_SOURCE_DIR}/rrtmgp/data-loading-examples)
#
# All-sky tests with clouds/by-band optics and aerosols
#
add_executable(rrtmgp_allsky rrtmgp_allsky.F90)
target_link_libraries(
  rrtmgp_allsky
  PRIVATE rrtmgp-load-data
          testing-utils
          rte-rrtmgp::rte
          rte-rrtmgp::rrtmgp
          rte-rrtmgp::ssm
)
# ##############################################################################
#
# RRTMGP tests and checks
#
add_test(
  NAME run_allsky_lw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw-bnd.nc
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
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw-bnd.nc
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

#
# All-sky tests with clouds/by-band optics only (no aerosols)
#
add_test(
  NAME run_allsky_no_aerosols_lw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw-bnd.nc
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
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw-bnd.nc
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

#
# All-sky tests with clouds/g-point optics only (no aerosols)
#
add_test(
  NAME run_allsky_gptclds_lw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols-gpt-clouds.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
    ${RRTMGP_DATA}/rrtmgp-clouds-lw-g256.nc
)
set_tests_properties(
  run_allsky_gptclds_lw
  PROPERTIES FIXTURES_REQUIRED
             fetch_rrtmgp_data
             FIXTURES_SETUP
             run_allsky_gpt_clouds
)

add_test(
  NAME run_allsky_gptclds_sw
  COMMAND
    rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw-no-aerosols-gpt-clouds.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
    ${RRTMGP_DATA}/rrtmgp-clouds-sw-g224.nc
)
set_tests_properties(
  run_allsky_gptclds_sw
  PROPERTIES FIXTURES_REQUIRED
             fetch_rrtmgp_data
             FIXTURES_SETUP
             run_allsky_gpt_clouds
)

add_test(
  NAME check_allsky_gptclds_lw_sw
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
    --ref_dir ${RRTMGP_DATA}/examples/all-sky/reference --tst_dir
    ${CMAKE_CURRENT_BINARY_DIR} --variables lw_flux_up lw_flux_dn sw_flux_up
    sw_flux_dn sw_flux_dir --file_names
    rrtmgp-allsky-lw-no-aerosols-gpt-clouds.nc
    rrtmgp-allsky-sw-no-aerosols-gpt-clouds.nc --failure_threshold
    ${FAILURE_THRESHOLD}
)
set_tests_properties(
  check_allsky_gptclds_lw_sw
  PROPERTIES FIXTURES_REQUIRED "fetch_rrtmgp_data;run_allsky_gpt_clouds"
)
# ##############################################################################
#
# SSM tests and checks
#
add_test(NAME link_ssm_allsky COMMAND ln rrtmgp_allsky ssm_allsky)
set_tests_properties(
  link_ssm_allsky PROPERTIES FIXTURES_SETUP link_allsky_executables
)

add_test(NAME run_ssm_allsky_lw COMMAND ssm_allsky 24 72 1 rrtmgp-allsky-lw.nc)
set_tests_properties(
  run_ssm_allsky_lw
  PROPERTIES FIXTURES_REQUIRED
             link_allsky_executables
             FIXTURES_SETUP
             run_rrtmgp_allsky
)

add_test(
  NAME run_ssm_allsky_sw COMMAND ssm_allsky 24 72 1 rrtmgp-allsky-sw.nc 5760
)
set_tests_properties(
  run_ssm_allsky_sw
  PROPERTIES FIXTURES_REQUIRED
             link_allsky_executables
             FIXTURES_SETUP
             run_rrtmgp_allsky
)
