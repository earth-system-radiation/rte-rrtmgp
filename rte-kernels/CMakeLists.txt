message(STATUS "Configuring lib rtekernels for '${KERNEL_MODE}' usacase")

# List of default files to compile.
# If KERNEL_MODE is accel or extern, then we try to compile
# accel versions or extern.
# If accel/extern versions of the files are not available, we fall back to the default
set(RTE_FORTRAN_KERNELS 
    mo_fluxes_broadband_kernels.F90
    mo_optical_props_kernels.F90
    mo_rte_solver_kernels.F90
    mo_rte_util_array.F90
)

# Build a list of files to compile based on input file names list
# For each file in the FILENAMES input array, we check if
# the path ${CMAKE_CURRENT_SOURCE_DIR}/${KERNEL_MODE}/{$FILENAME}
# is a valid path and if it is, we append the file to the output list
# If not, we append the default file instead.
function(find_files_in_folders OUTPUT_VAR FILENAMES KERNEL_MODE)
    set(RESULT_LIST) # Initialize an empty list to store results

    foreach(FILENAME IN LISTS FILENAMES)
        set(FULL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${KERNEL_MODE}/${FILENAME}")
        if(EXISTS "${FULL_PATH}")
            list(APPEND RESULT_LIST "${FULL_PATH}")
            message(STATUS "Using ${FULL_PATH}")
            continue()
        endif()
        set(FULL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
        list(APPEND RESULT_LIST "${FULL_PATH}")
        message(STATUS "Using ${FULL_PATH}")
    endforeach()

    # Return the result list to the caller
    set(${OUTPUT_VAR} "${RESULT_LIST}" PARENT_SCOPE)
endfunction()


# Variable to hold the output list
if (KERNEL_MODE STREQUAL "accel")
    set(MODE "accel")
elseif(KERNEL_MODE STREQUAL "extern")
    set(MODE "api")
elseif(KERNEL_MODE STREQUAL "default")
    set(MODE ".")
else()
    message(FATAL_ERROR "Invalid KERNEL_MODE ${KERNEL_MODE}")
endif()

find_files_in_folders(SOURCES "${RTE_FORTRAN_KERNELS}" "${MODE}")

add_library(rtekernels STATIC 
    ${SOURCES}
    ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_kind.F90
)