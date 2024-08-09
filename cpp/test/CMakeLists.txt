cmake_minimum_required(VERSION 3.0)
project(standalone C Fortran CXX)

enable_testing()

if ("${YAKL_ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

if (RRTMGP_ENABLE_YAKL)
  add_definitions("-DRRTMGP_ENABLE_YAKL")
endif()

if (RRTMGP_ENABLE_KOKKOS)
  add_definitions("-DRRTMGP_ENABLE_KOKKOS")
endif()

set (CMAKE_CXX_STANDARD 20)

# Build RRTMGP, YAKL, and the standalone allsky code (which depends on both)
# NOTE: ideally, we would build RRTMGP as a library and link allsky against it,
# but that does not seem to be working for GPU builds on summit right now for
# some reason (leads to segfaults). So instead, the allsky build will just pull
# in all the RRTMGP source and build and link all at once.
if (RRTMGP_ENABLE_YAKL)
  add_subdirectory(${YAKL_HOME} ./yakl)
endif()
if (RRTMGP_ENABLE_KOKKOS)
  add_subdirectory(${Kokkos_DIR} ./kokkos)
endif()
add_subdirectory(../ ./rrtmgp)
add_subdirectory(../examples/all-sky ./allsky)
add_subdirectory(./allsky_fortran ./allsky_fortran)
add_subdirectory(./allsky_fortran_openacc ./allsky_fortran_openacc)

# Add unit test script commands
add_test(NAME allsky_lw COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/build/test_lw.sh)
add_test(NAME allsky_sw COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/build/test_sw.sh)
add_test(NAME allsky_lw_r COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/build/test_lw_r.sh)
add_test(NAME allsky_sw_r COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/build/test_sw_r.sh)
