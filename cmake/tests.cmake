message(STATUS "[Tests enabled]")
message(STATUS "Configuring tests and examples")

find_package(Python REQUIRED COMPONENTS Interpreter)
find_package(NetCDF_Fortran REQUIRED)

# Define failure threshold for python tests to use
set(FAILURE_THRESHOLD
    7.e-4
    CACHE STRING "Default failure threshold for tests"
)

message(STATUS "CTest failure threshold: ${FAILURE_THRESHOLD}")

# Function to define a custom test for CTest Arguments: TEST_NAME - The name of
# the test in CTest ARGN      - The libraries that are going to be linked
function(build_test TEST_NAME)
  set(TEST_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.F90)
  set(TEST_LIBS ${ARGN})

  add_executable(${TEST_NAME} ${TEST_SOURCE})
  target_link_libraries(${TEST_NAME} PRIVATE ${TEST_LIBS})

  message(STATUS "Configuring test ${TEST_NAME}")
endfunction()

# Function to define a custom test for CTest Arguments: TEST_NAME    - The name
# of the test in CTest COMMAND_NAME - The command to execute the test
function(add_custom_test TEST_NAME COMMAND_NAME)
  # Add the test to CTest with a specific command and working directory
  add_test(
    NAME ${TEST_NAME}
    COMMAND ${COMMAND_NAME} ${ARGN}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  message(STATUS "Added test ${TEST_NAME}")
endfunction()

set(TEST_ATMOSPHERES
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
)
