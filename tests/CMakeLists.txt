SET(TEST_DIR ${CMAKE_SOURCE_DIR}/tests)
SET(RRTMGP_DATA ${CMAKE_SOURCE_DIR}/rrtmgp-data)

SET(LIBS rrtmgp rte netcdff)

# SET data directory for input files
SET(TEST_INPUT_FILES
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
)

# Additional sources for all tests
SET(ADDITIONAL_SOURCES
    ${TEST_DIR}/mo_testing_io.F90
    ${TEST_DIR}/mo_testing_utils.F90
    ${TEST_DIR}/mo_gas_optics_defs_rrtmgp.F90

    ${CMAKE_SOURCE_DIR}/examples/mo_simple_netcdf.F90
    ${CMAKE_SOURCE_DIR}/examples/rfmip-clear-sky/mo_rfmip_io.F90

    ${CMAKE_SOURCE_DIR}/extensions/mo_heating_rates.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_compute_bc.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_rrtmgp_clr_all_sky.F90
    ${CMAKE_SOURCE_DIR}/extensions/solar_variability/mo_solar_variability.F90
)

add_library(sources_objects OBJECT ${ADDITIONAL_SOURCES})
set_target_properties(sources_objects PROPERTIES Fortran_MODULE_DIRECTORY
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Check if netcdf.mod exists in /usr/local/include, then copy it
SET(NETCDF_DIR /usr/include)
if(EXISTS "${NETCDF_DIR}/netcdf.mod")
    file(COPY "${NETCDF_DIR}/netcdf.mod" DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY})
else()
    message(FATAL_ERROR "netcdf.mod not found in ${NETCDF_DIR}")
endif()

# Define each test
# foreach(TEST rte_optic_prop_unit_tests rte_lw_solver_unit_tests rte_sw_solver_unit_tests check_equivalence check_variants)
#     add_executable(${TEST} ${TEST_DIR}/${TEST}.F90 $<TARGET_OBJECTS:sources_objects>)
#     target_link_libraries(${TEST} PRIVATE ${LIBS})
#     message("-- Building test ${TEST}.so")
# endforeach()

add_executable(test_zenith_angle_spherical_correction
    ${CMAKE_SOURCE_DIR}/extensions/mo_zenith_angle_spherical_correction.F90
    $<TARGET_OBJECTS:sources_objects>
)
target_link_libraries(
    test_zenith_angle_spherical_correction PRIVATE ${LIBS})
message("-- Building test test_zenith_angle_spherical_correction.so")

# Copy the necessary input file for the tests
# add_custom_command(
#     TARGET test_zenith_angle_spherical_correction PRE_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy ${TEST_INPUT_FILES} ${CMAKE_CURRENT_BINARY_DIR}/test_atmospheres.nc
# )

# # Custom target to run all tests
# add_custom_target(run_tests
#     COMMAND ${CMAKE_COMMAND} -E copy ${TEST_INPUT_FILES} ${CMAKE_CURRENT_BINARY_DIR}/test_atmospheres.nc
#     COMMAND ${RUN_CMD} bash all_tests.sh
#     DEPENDS check_equivalence check_variants test_zenith_angle_spherical_correction
#     COMMENT "Running all tests"
# )
