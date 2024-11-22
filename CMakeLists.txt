cmake_minimum_required(VERSION 3.4...3.18)

project(rte-rrtmgp Fortran)

# Compiler settings.

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
    set(INCLUDE_DIRS "/usr/local/include")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(USER_Fortran_FLAGS "-Wall -fopenmp -fPIC")
    set(USER_Fortran_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")
    set(USER_Fortran_FLAGS_DEBUG "-g -O0")
    set(INCLUDE_DIRS "/usr/local/include")
endif()

# find_package(Python REQUIRED COMPONENTS Interpreter Development)

# Set the default build type to RELEASE.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release." FORCE)
else()
    string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
    set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
        "Choose the type of build, options are: None Debug Release." FORCE)
endif()

if (NOT HAS_CACHE)
    set(CMAKE_Fortran_FLAGS ${USER_Fortran_FLAGS} CACHE STRING
        "Flags used by the Fortran-compiler during all build types." FORCE)
    set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} CACHE STRING
        "Flags used by the Fortran-compiler during release build types." FORCE)
    set(CMAKE_Fortran_FLAGS_DEBUG ${USER_Fortran_FLAGS_DEBUG} CACHE STRING
        "Flags used by the Fortran-compiler during debug build types." FORCE)

    message(STATUS "Build Type: " ${CMAKE_BUILD_TYPE})
    set(HAS_CACHE TRUE CACHE BOOL "CMakeCache.txt created." FORCE)

    # Make sure that ccmake only contains build type.
    mark_as_advanced(HAS_CACHE)
endif()

add_subdirectory(rte-kernels)
add_subdirectory(rte-frontend)
