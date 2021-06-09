cmake_minimum_required(VERSION 3.0)
project(standalone)

# NOTE: Ideally, we would build rrtmgp as a library before building the
# the standalone code. However, this does not seem to be linking properly for
# GPU builds at the moment, so pull in all the RRTMGP source instead
set(RRTMGP_SRC
  ../../rrtmgp/kernels/mo_gas_optics_kernels.cpp
  ../../rrtmgp/mo_rrtmgp_constants.cpp
  ../../rrtmgp/mo_rrtmgp_util_reorder.cpp
  ../../rte/expand_and_transpose.cpp
  ../../rte/kernels/mo_fluxes_broadband_kernels.cpp
  ../../rte/kernels/mo_optical_props_kernels.cpp
  ../../rte/kernels/mo_rte_solver_kernels.cpp
)
set(CPP_SRC 
    ../mo_load_coefficients.cpp
    ../../extensions/fluxes_byband/mo_fluxes_byband_kernels.cpp
    mo_garand_atmos_io.cpp 
    mo_load_cloud_coefficients.cpp
    rrtmgp_allsky.cpp
    ${RRTMGP_SRC}
)
add_executable(allsky ${CPP_SRC})

target_link_libraries(allsky yakl ${NCFLAGS})
target_compile_features(allsky PUBLIC cxx_std_14)
target_include_directories(allsky PUBLIC ../)
target_include_directories(allsky PUBLIC ../../)
target_include_directories(allsky PUBLIC ../../extensions/cloud_optics)
target_include_directories(allsky PUBLIC ../../extensions/fluxes_byband)
target_include_directories(allsky PUBLIC ../../rte)
target_include_directories(allsky PUBLIC ../../rte/kernels)
target_include_directories(allsky PUBLIC ../../rrtmgp)
target_include_directories(allsky PUBLIC ../../rrtmgp/kernels)
target_include_directories(allsky PUBLIC ${NCINCLUDE})
set_source_files_properties( ${CPP_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS} ${CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(${CPP_SRC} PROPERTIES LANGUAGE CUDA)
endif()
