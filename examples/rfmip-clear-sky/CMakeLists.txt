set(TEST_INPUTS
    ${TEST_ATMOSPHERES}
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
)

# Copy required test inputs to the build directory
add_custom_target(
  rfmip_clear_sky_test_inputs ALL
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different ${TEST_INPUTS}
    ${CMAKE_CURRENT_BINARY_DIR}/
  COMMENT "Copying required .nc test files to build directory"
)

add_library(rfmip_clear_utils STATIC mo_rfmip_io.F90)
target_link_libraries(rfmip_clear_utils PRIVATE examples_utils)

foreach(
  test_executable IN
  ITEMS # cmake-format: sort
        rrtmgp_rfmip_lw
        rrtmgp_rfmip_sw
)
  add_executable(${test_executable} ${test_executable}.F90)
  target_link_libraries(${test_executable} rfmip_clear_utils examples_utils)
  add_dependencies(${test_executable} rfmip_clear_sky_test_inputs)
endforeach()

add_test(
  NAME rfmip_test_lw
  COMMAND
    rrtmgp_rfmip_lw 8
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
)

add_test(
  NAME rfmip_test_sw
  COMMAND
    rrtmgp_rfmip_sw 8
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
)

add_test(
  NAME check_rfmip_lw_sw
  COMMAND
    ${CMAKE_COMMAND} -E env ${Python_EXECUTABLE}
    ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py --ref_dir
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/reference --tst_dir
    ${CMAKE_CURRENT_BINARY_DIR} --variables rld rlu rsd rsu --file_names
    rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc --failure_threshold
    ${FAILURE_THRESHOLD}
)
