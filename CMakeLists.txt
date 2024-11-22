cmake_minimum_required(VERSION 3.4...3.18)
project(rte-rrtmgp LANGUAGES Fortran)

set(LINUX     FALSE)
set(APPLE_ARM FALSE)

if (UNIX AND NOT APPLE)
    set(LINUX TRUE)
elseif (APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    set(APPLE_ARM TRUE)
endif()

set(PREFERRED_FC_COMPILERS "gfortran" "ifort" "gfortran-10" "gfortran-11" "gfortran-12" "f77")

set(RTE_KERNEL_MODE "default" CACHE STRING "Select the kernel mode: default, accel, or extern")
set_property(CACHE RTE_KERNEL_MODE PROPERTY STRINGS "default" "accel" "extern")

set(PRECISION "double" CACHE STRING "Select the kernel mode: double or single")
set_property(CACHE PRECISION PROPERTY STRINGS "double" "single")

set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "Select the Fortran compiler: ${PREFERRED_FC_COMPILERS}")
set_property(CACHE CMAKE_Fortran_COMPILER PROPERTY STRINGS ${PREFERRED_FC_COMPILERS})

string(TOUPPER ${PRECISION} TMP_PRECISION)
add_definitions("-DUSE_${TMP_PRECISION}_PRECISION")

find_program(FORTRAN_COMPILER_FOUND ${CMAKE_Fortran_COMPILER})
if(NOT FORTRAN_COMPILER_FOUND)
    message(FATAL_ERROR "No suitable Fortran compiler found. Please install one of the following: ${PREFERRED_FC_COMPILERS}")
endif()

set(RTE_DIR           ${CMAKE_SOURCE_DIR}/rte-frontend)
set(RTE_KERNEL_DIR    ${CMAKE_SOURCE_DIR}/rte-kernels)
set(RRTMGP_DIR        ${CMAKE_SOURCE_DIR}/rrtmgp-frontend)
set(RRTMGP_KERNEL_DIR ${CMAKE_SOURCE_DIR}/rrtmgp-kernels)
set(TESTS_DIR         ${CMAKE_SOURCE_DIR}/tests)
set(RRTMGP_DATA       ${CMAKE_SOURCE_DIR}/rrtmgp-data)

include(FetchContent)

FetchContent_Declare(
    rrtmgp-data
    GIT_REPOSITORY https://github.com/earth-system-radiation/rrtmgp-data.git
    GIT_TAG origin/main
    SOURCE_DIR ${RRTMGP_DATA}
)
FetchContent_GetProperties(rrtmgp-data)

if (NOT rrtmgp-data_POPULATED)
    FetchContent_Populate(rrtmgp-data)
endif()

message(STATUS "RTE_KERNEL_MODE: ${RTE_KERNEL_MODE}")
message(STATUS "PRECISION: ${PRECISION} (double = 64-bits, single = 32-bits float)")
message(STATUS "FORTRAN COMPILER: ${CMAKE_Fortran_COMPILER}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(RTE_KIND_SRC ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_kind.F90)
add_library(rtef_kind STATIC ${RTE_KIND_SRC})

message(STATUS "Building lib rtef_kind.a")

add_subdirectory(${RTE_KERNEL_DIR})    # Build lib: rtekernels.a
add_subdirectory(${RTE_DIR})           # Build lib: rte.a
add_subdirectory(${RRTMGP_KERNEL_DIR}) # Build lib: rrtmgpkernels.a
add_subdirectory(${RRTMGP_DIR})        # Build lib: rrtmgp.a

# Enable testing with CTest
enable_testing()

add_subdirectory(${TESTS_DIR})         # Test & Check
