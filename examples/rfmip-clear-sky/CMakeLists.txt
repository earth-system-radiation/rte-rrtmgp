add_library(rfmip_clear_utils STATIC mo_rfmip_io.F90)

target_link_libraries(rfmip_clear_utils PUBLIC testing-utils optics-utils)

foreach(
  test_executable IN
  ITEMS # cmake-format: sort
        rrtmgp_rfmip_lw
        rrtmgp_rfmip_sw
)
  add_executable(${test_executable} ${test_executable}.F90)
  target_link_libraries(${test_executable} PRIVATE rfmip_clear_utils)
endforeach()

set(inoutputs
    rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
    rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
)

list(
  TRANSFORM inoutputs
  PREPEND ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/
          OUTPUT_VARIABLE inputs
)

# The tests write to the input files, therefore we copy them:
add_test(
  NAME copy_rrtmgp_rfmip_inputs
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different ${inputs}
    ${CMAKE_CURRENT_BINARY_DIR}/
)
set_tests_properties(
  copy_rrtmgp_rfmip_inputs
  PROPERTIES FIXTURES_REQUIRED
             fetch_rrtmgp_data
             FIXTURES_SETUP
             copy_rrtmgp_rfmip_inputs
)

add_test(
  NAME run_rrtmgp_rfmip_lw
  COMMAND
    rrtmgp_rfmip_lw 8
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
)
set_tests_properties(
  run_rrtmgp_rfmip_lw
  PROPERTIES FIXTURES_REQUIRED
             "fetch_rrtmgp_data;copy_rrtmgp_rfmip_inputs"
             FIXTURES_SETUP
             run_rrtmgp_rfmip
)

add_test(
  NAME run_rrtmgp_rfmip_sw
  COMMAND
    rrtmgp_rfmip_sw 8
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
)
set_tests_properties(
  run_rrtmgp_rfmip_sw
  PROPERTIES FIXTURES_REQUIRED
             "fetch_rrtmgp_data;copy_rrtmgp_rfmip_inputs"
             FIXTURES_SETUP
             run_rrtmgp_rfmip
)

add_test(
  NAME check_rfmip_lw_sw
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/examples/compare-to-reference.py
    --ref_dir ${RRTMGP_DATA}/examples/rfmip-clear-sky/reference --tst_dir
    ${CMAKE_CURRENT_BINARY_DIR} --variables rld rlu rsd rsu --file_names
    ${inoutputs} --failure_threshold ${FAILURE_THRESHOLD}
)
set_tests_properties(
  check_rfmip_lw_sw
  PROPERTIES FIXTURES_REQUIRED "fetch_rrtmgp_data;run_rrtmgp_rfmip"
)
# ##############################################################################
#
# Simple spectral model tests
#
if(BUILD_SSM)
  target_link_libraries(rfmip_clear_utils PUBLIC testing-utils optics-utils-ssm)

  foreach(
    test_executable IN
    ITEMS # cmake-format: sort
          rfmip_lw
          rfmip_sw
  )
    add_executable(ssm_${test_executable} rrtmgp_${test_executable}.F90)
    target_link_libraries(ssm_${test_executable} PRIVATE rfmip_clear_utils)
  endforeach()
  #
  # Copy and rename files
  #
  set(RENAMED_OUTPUTS)

  foreach(fn IN LISTS inoutputs)
    # Rename rule: p1f1 -> p1f2
    string(REGEX REPLACE "p1f1" "p1f2" new_name "${fn}")

    set(src "${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/${fn}")
    set(dst "${CMAKE_CURRENT_BINARY_DIR}/${new_name}")

    add_custom_command(
      OUTPUT "${dst}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "${src}" "${dst}"
      DEPENDS "${src}"
      COMMENT "Copying ${fn} → ${new_name}"
    )

    list(APPEND RENAMED_OUTPUTS "${dst}")
  endforeach()

  add_custom_target(
    copy_and_rename_rfmip_ssm_files ALL
    DEPENDS ${RENAMED_OUTPUTS}
    COMMENT "Copy and rename RFMIP -> SSM files"
  )

  add_test(
    NAME run_ssm_rfmip_lw
    COMMAND
      ssm_rfmip_lw 8
      ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
      ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc 2
  )
  set_tests_properties(
    run_ssm_rfmip_lw
    PROPERTIES FIXTURES_REQUIRED
               "fetch_rrtmgp_data;copy_rrtmgp_rfmip_inputs"
               FIXTURES_SETUP
               run_ssm_rfmip
               DEPENDS
               copy_and_rename_rfmip_ssm_files
  )

  add_test(
    NAME run_ssm_rfmip_sw
    COMMAND
      ssm_rfmip_sw 8
      ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
      ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc 2
  )
  set_tests_properties(
    run_ssm_rfmip_sw
    PROPERTIES FIXTURES_REQUIRED
               "fetch_rrtmgp_data;copy_rrtmgp_rfmip_inputs"
               FIXTURES_SETUP
               run_ssm_rfmip
               DEPENDS
               copy_and_rename_rfmip_ssm_files
  )

endif()
