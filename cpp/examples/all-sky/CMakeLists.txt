cmake_minimum_required(VERSION 3.0)
project(standalone)

set(CPP_SRC 
    ../mo_load_coefficients.cpp
    mo_garand_atmos_io.cpp 
    mo_load_cloud_coefficients.cpp
    rrtmgp_allsky.cpp
)
add_executable(allsky ${CPP_SRC})

target_link_libraries(allsky yakl rrtmgp ${NCFLAGS})
target_compile_features(allsky PUBLIC cxx_std_14)

target_include_directories(allsky PUBLIC ${YAKL_HOME})
target_include_directories(allsky PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/yakl)
set_source_files_properties( ${CPP_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS} ${CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
  set_source_files_properties(${CPP_SRC} PROPERTIES LANGUAGE CUDA)
  target_include_directories(allsky PUBLIC ${YAKL_CUB_HOME})
endif()
