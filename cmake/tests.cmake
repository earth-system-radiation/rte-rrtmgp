# Define failure threshold for python tests to use
set(FAILURE_THRESHOLD
    7.e-4
    CACHE STRING "Default failure threshold for tests"
)

message(STATUS "CTest failure threshold: ${FAILURE_THRESHOLD}")

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
