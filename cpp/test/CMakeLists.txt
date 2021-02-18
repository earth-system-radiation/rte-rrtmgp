cmake_minimum_required(VERSION 3.0)
project(standalone)

enable_testing()

if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

# Build RRTMGP, YAKL, and the standalone allsky code (which depends on both)
# NOTE: ideally, we would build RRTMGP as a library and link allsky against it,
# but that does not seem to be working for GPU builds on summit right now for
# some reason (leads to segfaults). So instead, the allsky build will just pull
# in all the RRTMGP source and build and link all at once.
#add_subdirectory(../ ./rrtmgp)
add_subdirectory(${YAKL_HOME} ./yakl)
add_subdirectory(../examples/all-sky ./allsky)

# Add unit test script commands
add_test(NAME allsky_lw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_lw.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME allsky_sw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_sw.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME allsky_lw_r COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_lw_r.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME allsky_sw_r COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_sw_r.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
