SET(TEST_DIR ${CMAKE_SOURCE_DIR}/tests)
SET(CTEST_TEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Path to the test atmospheres file used for equivalence and variant tests
SET(TEST_ATMOSPHERES ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc)

# Function to build a test
# Arguments:
#   TEST_NAME - The name of the test executable to create
function(build_test TEST_NAME)
    set(TEST_SOURCE ${TEST_DIR}/${TEST_NAME}.F90)

    add_executable(${TEST_NAME} ${TEST_SOURCE} ${ARGN})
    target_link_libraries(${TEST_NAME} PRIVATE test_utils)

    message(STATUS "Building test ${TEST_NAME}.so")
endfunction()

# Function to define a custom test for CTest
# Arguments:
#   TEST_NAME    - The name of the test in CTest
#   COMMAND_NAME - The command to execute the test
function(add_custom_test TEST_NAME COMMAND_NAME)
    # Add the test to CTest with a specific command and working directory
    add_test(NAME ${TEST_NAME}
        COMMAND ${COMMAND_NAME} ${ARGN}
        WORKING_DIRECTORY ${CTEST_TEST_DIR}
    )

    message(STATUS "Added test ${TEST_NAME} with command ${COMMAND_NAME}")
endfunction()

# List of additional source files shared across tests
SET(SOURCES
    ${TEST_DIR}/mo_testing_io.F90
    ${TEST_DIR}/mo_testing_utils.F90
    ${TEST_DIR}/mo_rcemip_profiles.F90
    ${TEST_DIR}/mo_gas_optics_defs_rrtmgp.F90

    ${CMAKE_SOURCE_DIR}/examples/mo_simple_netcdf.F90
    ${CMAKE_SOURCE_DIR}/examples/rfmip-clear-sky/mo_rfmip_io.F90

    ${CMAKE_SOURCE_DIR}/extensions/mo_compute_bc.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_heating_rates.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_rrtmgp_clr_all_sky.F90
    ${CMAKE_SOURCE_DIR}/extensions/mo_zenith_angle_spherical_correction.F90
    ${CMAKE_SOURCE_DIR}/extensions/solar_variability/mo_solar_variability.F90
)

# Check if netcdf.mod exists in /usr/local/include, then copy it
# TODO: Find a better way to use netcdf.mod

find_package(NetCDF COMPONENTS C Fortran REQUIRED)
message(STATUS "=============================== NETCDF")
# List the available targets
message(STATUS "Available NetCDF targets:")
message(STATUS "  NetCDF::NetCDF -> ${NetCDF_LIBRARIES}")
message(STATUS "  NetCDF::NetCDF_Fortran -> ${NetCDF_Fortran_LIBRARIES}")
message(STATUS "=============================== NETCDF")

if (NOT WIN32)
    set(NETCDF_DIR /usr/include)
    # Check if netcdf.mod exists in the specified directory
    if(EXISTS "${NETCDF_DIR}/netcdf.mod")
        # Copy netcdf.mod to the Fortran module directory
        file(COPY "${NETCDF_DIR}/netcdf.mod" DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY})
    else()
        message(SEND_ERROR "netcdf.mod not found in ${NETCDF_DIR}.")

        if(LINUX)
            message(FATAL_ERROR "On Linux, you can install it using: `sudo apt-get install libnetcdff-dev`")
        elseif(APPLE_ARM)
            message(FATAL_ERROR "On macOS ARM (Apple Silicon), you can install it using: `brew install netcdf`")
        else()
            message(FATAL_ERROR "Please ensure NetCDF Fortran is installed and netcdf.mod is in ${NETCDF_DIR}.")
        endif()
    endif()
endif()

add_library(test_utils STATIC ${SOURCES})

# Link necessary libraries to the test utility library
target_link_libraries(test_utils PUBLIC
    rte
    rrtmgp
    netCDF::netcdf
)

# Set the directory for Fortran module files
set_target_properties(test_utils PROPERTIES Fortran_MODULE_DIRECTORY
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Build tests
build_test(rte_optic_prop_unit_tests)
build_test(rte_lw_solver_unit_tests)
build_test(rte_sw_solver_unit_tests)

build_test(check_equivalence)
build_test(check_variants)

# No tests are provided at the moment
build_test(test_zenith_angle_spherical_correction)

# Add unit tests
add_custom_test(rte_optic_prop_unit_tests rte_optic_prop_unit_tests)
add_custom_test(rte_lw_solver_unit_tests  rte_lw_solver_unit_tests)
add_custom_test(rte_sw_solver_unit_tests  rte_sw_solver_unit_tests)

# Add equivalence tests
add_custom_test(check_equivalence_lw_g256 check_equivalence ${TEST_ATMOSPHERES} ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc)
add_custom_test(check_equivalence_lw_g128 check_equivalence ${TEST_ATMOSPHERES} ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc)
add_custom_test(check_equivalence_sw_g224 check_equivalence ${TEST_ATMOSPHERES} ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc)
add_custom_test(check_equivalence_sw_g112 check_equivalence ${TEST_ATMOSPHERES} ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc)

# Add variant tests
add_custom_test(check_variants_lw check_variants ${TEST_ATMOSPHERES} ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc)
add_custom_test(check_variants_sw check_variants ${TEST_ATMOSPHERES} ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc)
