cmake_minimum_required(VERSION 3.0)
project(standalone)

enable_testing()

if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

file(GLOB CPP_SRC ../rrtmgp/*.cpp
                  ../rrtmgp/kernels/*.cpp
                  ../rte/*.cpp
                  ../rte/kernels/*.cpp
                  ../extensions/cloud_optics/mo_cloud_optics.cpp
                  ../examples/*.cpp )

add_subdirectory(${YAKL_HOME} ./yakl)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/yakl)

add_subdirectory(allsky)

