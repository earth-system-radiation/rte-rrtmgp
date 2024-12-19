message(STATUS "Configuring lib examples_utils")

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_garand_atmos_io.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_load_coefficients.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/mo_simple_netcdf.F90
)

add_library(examples_utils STATIC ${SOURCES})

target_include_directories(examples_utils PUBLIC ${NetCDF_Fortran_INCLUDE_DIR})

target_link_libraries(
  examples_utils PUBLIC rte rrtmgp ${NetCDF_Fortran_LIBRARY}
)

add_subdirectory(all-sky)
add_subdirectory(rfmip-clear-sky)
