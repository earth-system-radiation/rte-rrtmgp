cmake_minimum_required(VERSION 3.4...3.18)
project(rte-rrtmgp LANGUAGES Fortran)

set(RTE_KERNEL_MODE "default" CACHE STRING "Select the kernel mode: default, accel, or extern")
set_property(CACHE RTE_KERNEL_MODE PROPERTY STRINGS "default" "accel" "extern")

set(PRECISION "double" CACHE STRING "Select the kernel mode: double or single")
set_property(CACHE PRECISION PROPERTY STRINGS "double" "single")

string(TOUPPER ${PRECISION} TMP_PRECISION)
add_definitions("-DUSE_${TMP_PRECISION}_PRECISION")

message(STATUS "RTE_KERNEL_MODE: ${RTE_KERNEL_MODE}")
message(STATUS "PRECISION: ${PRECISION} (double = 64-bits, single = 32-bits float)")

# Fortran compiler selection
set(PREFERRED_FC_COMPILERS gfortran ifort gfortran-10 gfortran-11 gfortran-12 f77)

foreach(compiler IN LISTS PREFERRED_FC_COMPILERS)
    find_program(FOUND_COMPILER NAMES ${compiler})
    if(FOUND_COMPILER)
        set(CMAKE_Fortran_COMPILER ${FOUND_COMPILER})
        message(STATUS "Using Fortran compiler: ${FOUND_COMPILER}")
        break()
    endif()
endforeach()

if(NOT CMAKE_Fortran_COMPILER)
    message(FATAL_ERROR "No suitable Fortran compiler found")
endif()

SET(RTE_DIR           ${CMAKE_SOURCE_DIR}/rte-frontend)
SET(RTE_KERNEL_DIR    ${CMAKE_SOURCE_DIR}/rte-kernels)
SET(RRTMGP_DIR        ${CMAKE_SOURCE_DIR}/rrtmgp-frontend)
SET(RRTMGP_KERNEL_DIR ${CMAKE_SOURCE_DIR}/rrtmgp-kernels)
SET(TESTS_DIR         ${CMAKE_SOURCE_DIR}/tests)
SET(RRTMGP_DATA       ${CMAKE_SOURCE_DIR}/rrtmgp-data)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

set(RTE_KIND_SRC ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_kind.F90)
add_library(rtef_kind STATIC ${RTE_KIND_SRC})

message(STATUS "Building lib rtef_kind.a")

add_subdirectory(${RTE_KERNEL_DIR})    # Build lib: rtekernels.a
add_subdirectory(${RTE_DIR})           # Build lib: rte.a
add_subdirectory(${RRTMGP_KERNEL_DIR}) # Build lib: rrtmgpkernels.a
add_subdirectory(${RRTMGP_DIR})        # Build lib: rrtmgp.a
# add_subdirectory(${TESTS_DIR})         # Test & Check
