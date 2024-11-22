set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_garand_atmos_io.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_load_coefficients.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_simple_netcdf.F90
)

add_library(examples_utils STATIC ${SOURCES})

target_include_directories(examples_utils PUBLIC
    ${NETCDF_Fortran_INCLUDE_DIR}
)

target_link_libraries(examples_utils PUBLIC
    rte
    rrtmgp
    ${NETCDF_Fortran_LIBRARY}
)

message(STATUS "Building lib examples_utils.a")

find_package(Python REQUIRED COMPONENTS Interpreter)

# Define failure threshold as an environment variable
set(FAILURE_THRESHOLD 7.e-4 CACHE STRING "Default failure threshold for tests")

# TODO: Fix all-sky
add_subdirectory(all-sky)
add_subdirectory(rfmip-clear-sky)
