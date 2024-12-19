set(TEST_INPUTS
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
)

# Copy required test inputs to the build directory
add_custom_target(
  tests_test_inputs ALL
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different ${TEST_INPUTS}
    ${CMAKE_CURRENT_BINARY_DIR}/
  COMMENT "Copying required .nc test files to build directory"
)

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

add_dependencies(check_equivalence tests_test_inputs)
add_dependencies(check_variants tests_test_inputs)

add_test(NAME rte_optic_prop_unit_tests COMMAND rte_optic_prop_unit_tests)
add_test(NAME rte_lw_solver_unit_tests COMMAND rte_lw_solver_unit_tests)
add_test(NAME rte_sw_solver_unit_tests COMMAND rte_sw_solver_unit_tests)

add_test(
  NAME check_equivalence_lw_g256
  COMMAND
    check_equivalence
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
)

add_test(
  NAME check_equivalence_lw_g128
  COMMAND
    check_equivalence
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)

add_test(
  NAME check_equivalence_sw_g224
  COMMAND
    check_equivalence
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
)

add_test(
  NAME check_equivalence_sw_g112
  COMMAND
    check_equivalence
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)

add_test(
  NAME check_variants_lw
  COMMAND
    check_variants
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)

add_test(
  NAME check_variants_sw
  COMMAND
    check_variants
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)

add_custom_target(
  validation-plots
  COMMAND
    ${CMAKE_COMMAND} -E env ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/validation-plots.py --input_file
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  COMMENT "Generating validation plots using validation-plots.py"
)
