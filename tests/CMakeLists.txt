set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

set(extensions_source_dir ${PROJECT_SOURCE_DIR}/extensions)

add_library(
  test_utils STATIC # cmake-format: sort
  ${extensions_source_dir}/mo_compute_bc.F90
  ${extensions_source_dir}/mo_heating_rates.F90
  ${extensions_source_dir}/mo_rrtmgp_clr_all_sky.F90
  ${extensions_source_dir}/mo_zenith_angle_spherical_correction.F90
  ${extensions_source_dir}/solar_variability/mo_solar_variability.F90
  mo_gas_optics_defs_rrtmgp.F90
  mo_rcemip_profiles.F90
  mo_testing_utils.F90
)

target_include_directories(
  test_utils
  PUBLIC
    $<BUILD_INTERFACE:$<$<COMPILE_LANGUAGE:Fortran>:${CMAKE_Fortran_MODULE_DIRECTORY}>>
    ${NetCDF_Fortran_INCLUDE_DIR}
)

target_link_libraries(test_utils PUBLIC rfmip_clear_utils)

foreach(
  test_executable IN
  ITEMS # cmake-format: sort
        check_equivalence
        check_variants
        rte_lw_solver_unit_tests
        rte_optic_prop_unit_tests
        rte_sw_solver_unit_tests
        test_zenith_angle_spherical_correction
)
  add_executable(${test_executable} ${test_executable}.F90)
  target_link_libraries(${test_executable} PRIVATE test_utils)
endforeach()

foreach(
  test_executable IN
  ITEMS # cmake-format: sort
        rte_lw_solver_unit_tests
        rte_optic_prop_unit_tests
        rte_sw_solver_unit_tests
)
  add_test(NAME ${test_executable} COMMAND ${test_executable})
endforeach()

foreach(g_value IN ITEMS 256 128)
  add_test(
    NAME check_equivalence_lw_g${g_value}
    COMMAND
      check_equivalence
      ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
      ${RRTMGP_DATA}/rrtmgp-gas-lw-g${g_value}.nc
  )
  set_tests_properties(
    check_equivalence_lw_g${g_value}
    PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data
  )
endforeach()

foreach(g_value IN ITEMS 224 112)
  add_test(
    NAME check_equivalence_sw_g${g_value}
    COMMAND
      check_equivalence
      ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
      ${RRTMGP_DATA}/rrtmgp-gas-sw-g${g_value}.nc
  )
  set_tests_properties(
    check_equivalence_sw_g${g_value}
    PROPERTIES FIXTURES_REQUIRED fetch_rrtmgp_data
  )
endforeach()

add_test(
  NAME check_variants_lw
  COMMAND
    check_variants
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    test_atmospheres.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)

add_test(
  NAME check_variants_sw
  COMMAND
    check_variants
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    test_atmospheres.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)
set_tests_properties(
  check_variants_sw
  PROPERTIES # The test writes to the same file as check_variants_lw, therefore
             # we cannot run them in parallel:
             DEPENDS
             check_variants_lw
)

add_custom_target(
  validation-plots
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/validation-plots.py
    --input_file multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  COMMENT "Generating validation plots"
)
