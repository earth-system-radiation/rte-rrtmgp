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
  mo_testing_io.F90
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

add_test(NAME rte_optic_prop_unit_tests COMMAND rte_optic_prop_unit_tests)
add_test(NAME rte_lw_solver_unit_tests COMMAND rte_lw_solver_unit_tests)
add_test(NAME rte_sw_solver_unit_tests COMMAND rte_sw_solver_unit_tests)

add_test(
  NAME check_equivalence_lw_g256
  COMMAND
    check_equivalence
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
)

add_test(
  NAME check_equivalence_lw_g128
  COMMAND
    check_equivalence
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)

add_test(
  NAME check_equivalence_sw_g224
  COMMAND
    check_equivalence
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
)

add_test(
  NAME check_equivalence_sw_g112
  COMMAND
    check_equivalence
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)

# The tests write to the input file, therefore we copy it:
add_test(
  NAME copy_check_variants_input
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${CMAKE_CURRENT_BINARY_DIR}/
)
set_tests_properties(
  copy_check_variants_input PROPERTIES FIXTURES_SETUP copy_check_variants_input
)

add_test(
  NAME check_variants_lw
  COMMAND
    check_variants
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)
set_tests_properties(
  check_variants_lw
  PROPERTIES FIXTURES_REQUIRED
             copy_check_variants_input
             FIXTURES_SETUP
             check_variants_lw
)

add_test(
  NAME check_variants_sw
  COMMAND
    check_variants
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)
set_tests_properties(
  check_variants_sw
  # The test writes to the same file as check_variants_lw, therefore we must not
  # run them in parallel:
  PROPERTIES FIXTURES_REQUIRED check_variants_lw
)

add_custom_target(
  validation-plots
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/validation-plots.py
    --input_file multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  COMMENT "Generating validation plots"
)
