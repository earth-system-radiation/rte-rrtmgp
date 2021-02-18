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
add_library(rrtmgp ${CXX_SRC})
target_link_libraries(rrtmgp yakl)

# Special compile flags for C++ source
set_source_files_properties(${CXX_SRC} PROPERTIES COMPILE_FLAGS "${CPPDEFS} ${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(${CXX_SRC} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${CXX_SRC} PROPERTIES COMPILE_FLAGS "${CUDA_FLAGS}")
  set_target_properties(rrtmgp PROPERTIES LINKER_LANGUAGE CUDA CUDA_SEPARABLE_COMPILATION OFF CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()
target_include_directories(rrtmgp PUBLIC ${YAKL_BIN})
target_include_directories(rrtmgp PUBLIC ${YAKL_HOME})
target_include_directories(rrtmgp PUBLIC ./)
target_include_directories(rrtmgp PUBLIC rte)
target_include_directories(rrtmgp PUBLIC rte/kernels)
target_include_directories(rrtmgp PUBLIC rrtmgp)
target_include_directories(rrtmgp PUBLIC rrtmgp/kernels)
