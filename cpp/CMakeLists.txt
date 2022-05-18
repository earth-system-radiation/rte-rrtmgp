###############################################################################
## This assumes you have already run add_subdirectory(${YAKL_HOME} ${YAKL_BIN})
## It also assumes you've enabled CXX and C as languages
###############################################################################

# Add RRTMGP source files.
set(CXX_SRC
  rrtmgp/kernels/mo_gas_optics_kernels.cpp
  rrtmgp/mo_rrtmgp_constants.cpp
  rrtmgp/mo_rrtmgp_util_reorder.cpp
  rte/expand_and_transpose.cpp
  rte/kernels/mo_fluxes_broadband_kernels.cpp
  rte/kernels/mo_optical_props_kernels.cpp
  rte/kernels/mo_rte_solver_kernels.cpp
)

# Add RRTMGP as a library and make sure we link against YAKL
add_library(rrtmgp STATIC ${CXX_SRC})
target_link_libraries(rrtmgp yakl)

# Special compile flags for C++ source
include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(rrtmgp)
message(STATUS "YAKL Flags: ${YAKL_COMPILER_FLAGS}")
target_include_directories(rrtmgp PUBLIC ${YAKL_HOME})
target_include_directories(rrtmgp PUBLIC ./)
target_include_directories(rrtmgp PUBLIC rte)
target_include_directories(rrtmgp PUBLIC rte/kernels)
target_include_directories(rrtmgp PUBLIC rrtmgp)
target_include_directories(rrtmgp PUBLIC rrtmgp/kernels)
