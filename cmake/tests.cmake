# Define failure threshold for python tests to use
set(FAILURE_THRESHOLD
    7.e-4
    CACHE STRING "Default failure threshold for tests"
)

message(STATUS "CTest failure threshold: ${FAILURE_THRESHOLD}")

set(TEST_ATMOSPHERES
    ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc
)
