set(TEST_INPUTS ${TEST_ATMOSPHERES})

# Copy required test inputs to the build directory
add_custom_target(
  tests_test_inputs ALL
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different ${TEST_INPUTS}
    ${CMAKE_CURRENT_BINARY_DIR}/
  COMMENT "Copying required .nc test files to build directory"
)

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
target_include_directories(test_utils PRIVATE ${NetCDF_Fortran_INCLUDE_DIR})
target_link_libraries(
  test_utils
  PRIVATE rte rrtmgp examples_utils rfmip_clear_utils ${NetCDF_Fortran_LIBRARY}
)

set_target_properties(
  test_utils
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

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
  target_link_libraries(${test_executable} test_utils)
endforeach()

add_dependencies(check_equivalence tests_test_inputs)
add_dependencies(check_variants tests_test_inputs)

# Add unit tests
add_custom_test(rte_optic_prop_unit_tests rte_optic_prop_unit_tests)
add_custom_test(rte_lw_solver_unit_tests rte_lw_solver_unit_tests)
add_custom_test(rte_sw_solver_unit_tests rte_sw_solver_unit_tests)

# Add equivalence tests
add_custom_test(
  check_equivalence_lw_g256
  check_equivalence
  multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
)
add_custom_test(
  check_equivalence_lw_g128
  check_equivalence
  multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)
add_custom_test(
  check_equivalence_sw_g224
  check_equivalence
  multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
)
add_custom_test(
  check_equivalence_sw_g112
  check_equivalence
  multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)

# Add variant tests
add_custom_test(
  check_variants_lw
  check_variants
  multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc
  ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
)
add_custom_test(
  check_variants_sw
  check_variants
  multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc
  ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc
)

# Add a custom target to generate validation plots
add_custom_target(
  validation-plots
  COMMAND
    ${CMAKE_COMMAND} -E env ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/validation-plots.py --input_file
    multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating validation plots using validation-plots.py"
)
