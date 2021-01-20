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
  examples/mo_load_coefficients.cpp
  # This is a hack; I should not have to build this here, but we need it for the standalone RRTMGP
  # unit tests, and I could not figure out how to get it to build without adding it to the RRTMGP
  # library. TODO: revisit this!
  examples/all-sky/mo_garand_atmos_io.cpp
  examples/all-sky/mo_load_cloud_coefficients.cpp
  extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
)

# Add RRTMGP as a library and make sure we link against YAKL
add_library(rrtmgp ${CXX_SRC})
target_link_libraries(rrtmgp yakl)

# Special compile flags for C++ source
set_source_files_properties(${CXX_SRC} PROPERTIES COMPILE_FLAGS "${CPPDEFS} ${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(${CXX_SRC} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${CXX_SRC} PROPERTIES COMPILE_FLAGS "${CUDA_FLAGS}")
  target_include_directories(rrtmgp PUBLIC ${YAKL_CUB_HOME})
endif()
target_include_directories(rrtmgp PUBLIC ${YAKL_BIN})
target_include_directories(rrtmgp PUBLIC ${YAKL_HOME})
target_include_directories(rrtmgp PUBLIC ./)
target_include_directories(rrtmgp PUBLIC rte)
target_include_directories(rrtmgp PUBLIC rte/kernels)
target_include_directories(rrtmgp PUBLIC rrtmgp)
target_include_directories(rrtmgp PUBLIC rrtmgp/kernels)
target_include_directories(rrtmgp PUBLIC examples)
# Again, this is a hack to get the standalone tests to build that need the mo_garand_atmos_io.h header
target_include_directories(rrtmgp PUBLIC examples/all-sky)
target_include_directories(rrtmgp PUBLIC extensions/fluxes_byband)
target_include_directories(rrtmgp PUBLIC extensions/cloud_optics)
