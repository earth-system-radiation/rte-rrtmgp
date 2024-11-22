cmake_minimum_required(VERSION 3.4...3.18)

project(rte-rrtmgp Fortran)

# Helper function to check if a value is in a list
function(check_valid_option variable_name value valid_options)
    list(FIND ${valid_options} ${value} index)
    if (index EQUAL -1)
        message(FATAL_ERROR
            "${variable_name} has an invalid value: '${value}'.\n"
            "Valid options are: ${${valid_options}}.")
    endif()
endfunction()

set(VALID_RTE_KERNELS "ACCEL" "EXTERN" "DEFAULT")
set(RTE_KERNELS "DEFAULT")

# Set the RTE_KERNELS value to upper case only
string(TOUPPER "${RTE_KERNELS}" RTE_KERNELS)

# Validate RTE_KERNELS
check_valid_option("RTE_KERNELS" "${RTE_KERNELS}" "VALID_RTE_KERNELS")

message(STATUS "RTE_KERNELS: ${RTE_KERNELS}")

SET(RTE_DIR           ${CMAKE_SOURCE_DIR}/rte-frontend)
SET(RTE_KERNEL_DIR    ${CMAKE_SOURCE_DIR}/rte-kernels)
SET(RRTMGP_DIR        ${CMAKE_SOURCE_DIR}/rrtmgp-frontend)
SET(RRTMGP_KERNEL_DIR ${CMAKE_SOURCE_DIR}/rrtmgp-kernels)
SET(TESTS_DIR         ${CMAKE_SOURCE_DIR}/tests)

SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Compiler settings
if (NOT USESP)
    set(USESP FALSE)
endif()

if (USESP)
    add_definitions("-DUSE_SINGLE_PRECISION")
    message(STATUS "Precision: Single (32-bits floats)")
else()
    message(STATUS "Precision: Double (64-bits floats)")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "AppleClang")
    set(USER_Fortran_FLAGS "-Wall -fPIC")
    set(USER_Fortran_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")
    set(USER_Fortran_FLAGS_DEBUG "-g -O0")
    set(USER_Fortran_INCLUDE_DIRS "/usr/local/include")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(USER_Fortran_FLAGS "-ffree-line-length-none -m64 -march=native -DRTE_USE_CBOOL")
    set(USER_Fortran_FLAGS_RELEASE "-std=f2008 -fbounds-check -fmodule-private -fimplicit-none -finit-real=nan -g")
    set(USER_Fortran_FLAGS_DEBUG "-std=f2008 -fbounds-check -fmodule-private -fimplicit-none -finit-real=nan -g")
    set(USER_Fortran_INCLUDE_DIRS "/usr/local/include")
endif()

# Set default build type to RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release." FORCE)
else()
    string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
    set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
        "Choose the type of build, options are: None Debug Release." FORCE)
endif()

# Set Fortran compiler flags if not already cached
if (NOT HAS_CACHE)
    set(CMAKE_Fortran_FLAGS ${USER_Fortran_FLAGS} CACHE STRING
        "Flags used by the Fortran-compiler during all build types." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} CACHE STRING
        "Flags used by the Fortran-compiler during release build types." FORCE)
    set(CMAKE_Fortran_FLAGS_DEBUG ${USER_Fortran_FLAGS_DEBUG} CACHE STRING
        "Flags used by the Fortran-compiler during debug build types." FORCE)

    message(STATUS "Build Type: " ${CMAKE_BUILD_TYPE})
    set(HAS_CACHE TRUE CACHE BOOL "CMakeCache.txt created." FORCE)
    mark_as_advanced(HAS_CACHE)
endif()

set(RTE_KIND_SRC ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_kind.F90)
add_library(rtef_kind STATIC ${RTE_KIND_SRC})

message(STATUS "Building lib rtef_kind.a")

add_subdirectory(${RTE_KERNEL_DIR})    # Build lib: rtekernels.a
add_subdirectory(${RTE_DIR})           # Build lib: rte.a
add_subdirectory(${RRTMGP_KERNEL_DIR}) # Build lib: rrtmgpkernels.a
add_subdirectory(${RRTMGP_DIR})        # Build lib: rrtmgp.a
# add_subdirectory(${TESTS_DIR})         # Test & Check
