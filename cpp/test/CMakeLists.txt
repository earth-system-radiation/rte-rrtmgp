cmake_minimum_required(VERSION 3.0)
project(standalone)

enable_testing()

if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

add_subdirectory(../ ./rrtmgp)
add_subdirectory(${YAKL_HOME} ./yakl)
add_subdirectory(../examples/all-sky ./allsky)

add_test(NAME allsky_lw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_lw.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME allsky_sw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_sw.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME allsky_lw_r COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_lw_r.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME allsky_sw_r COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_sw_r.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
