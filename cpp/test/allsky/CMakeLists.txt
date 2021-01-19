
set(SRCLOC ../../examples/all-sky/mo_garand_atmos_io.cpp 
           ../../examples/all-sky/mo_load_cloud_coefficients.cpp
           ../../examples/mo_load_coefficients.cpp
           ../../examples/all-sky/rrtmgp_allsky.cpp
           ${CPP_SRC})

add_executable(allsky ${SRCLOC})
target_link_libraries(allsky yakl ${NCFLAGS})
target_compile_features(allsky PUBLIC cxx_std_14)

set_source_files_properties( ${SRCLOC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS} ${CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${SRCLOC} PROPERTIES LANGUAGE CUDA)
endif()

include_directories(${YAKL_HOME})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../..)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../rte/kernels)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../rrtmgp/kernels)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../rte)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../extensions)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../extensions/cloud_optics)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../examples)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../rrtmgp)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/yakl)
include_directories(${YAKL_CUB_HOME})

add_test(NAME allsky_lw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_lw.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
add_test(NAME allsky_sw COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_sw.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
add_test(NAME allsky_lw_r COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_lw_r.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
add_test(NAME allsky_sw_r COMMAND ${CMAKE_CURRENT_BINARY_DIR}/../test_sw_r.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)
