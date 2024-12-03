set(TEST_INPUTS
    ${TEST_ATMOSPHERES}
)

# Copy required test inputs to the build directory
add_custom_target(tests_test_inputs ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${TEST_INPUTS} ${CMAKE_CURRENT_BINARY_DIR}/
    COMMENT "Copying required .nc test files to build directory"
)

# List of additional source files shared across tests
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_testing_io.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_testing_utils.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_rcemip_profiles.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_gas_optics_defs_rrtmgp.F90

    ${CMAKE_SOURCE_DIR}/extensions/mo_compute_bc.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_heating_rates.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_rrtmgp_clr_all_sky.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_zenith_angle_spherical_correction.F90
    ${CMAKE_SOURCE_DIR}/extensions/solar_variability/mo_solar_variability.F90
)

add_library(test_utils STATIC ${SOURCES})

target_include_directories(test_utils PRIVATE ${NETCDF_Fortran_INCLUDE_DIR})

target_link_libraries(test_utils PRIVATE
    rte
    rrtmgp
    examples_utils
    rfmip_clear_utils
    ${NETCDF_Fortran_LIBRARY}
)

# Set the directory for Fortran module files
set_target_properties(test_utils PROPERTIES Fortran_MODULE_DIRECTORY
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Build tests
build_test(rte_optic_prop_unit_tests test_utils)
build_test(rte_lw_solver_unit_tests  test_utils)
build_test(rte_sw_solver_unit_tests  test_utils)

build_test(check_equivalence test_utils)
build_test(check_variants    test_utils)

# No tests are provided at the moment
build_test(test_zenith_angle_spherical_correction test_utils)

add_dependencies(check_equivalence tests_test_inputs)
add_dependencies(check_variants tests_test_inputs)

# Add unit tests
add_custom_test(rte_optic_prop_unit_tests rte_optic_prop_unit_tests)
add_custom_test(rte_lw_solver_unit_tests  rte_lw_solver_unit_tests)
add_custom_test(rte_sw_solver_unit_tests  rte_sw_solver_unit_tests)

# Add equivalence tests
add_custom_test(check_equivalence_lw_g256 check_equivalence multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc)
add_custom_test(check_equivalence_lw_g128 check_equivalence multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc)
add_custom_test(check_equivalence_sw_g224 check_equivalence multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc)
add_custom_test(check_equivalence_sw_g112 check_equivalence multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc)

# Add variant tests
add_custom_test(check_variants_lw check_variants multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc)
add_custom_test(check_variants_sw check_variants multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc)

# Add a custom target to generate validation plots
add_custom_target(validation-plots
    COMMAND ${CMAKE_COMMAND} -E env ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/validation-plots.py
        --input_file multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating validation plots using validation-plots.py"
)
