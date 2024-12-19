message(STATUS "Configuring lib rtekernels for '${KERNEL_MODE}' usacase")

# List of default files to compile. If KERNEL_MODE is accel or extern, then we
# try to compile accel versions or extern. If accel/extern versions of the files
# are not available, we fall back to the default
set(RTE_FORTRAN_KERNELS
    mo_fluxes_broadband_kernels.F90
    mo_optical_props_kernels.F90
    mo_rte_solver_kernels.F90
    mo_rte_util_array.F90
)

find_files_in_folders(SOURCES "${RTE_FORTRAN_KERNELS}" "${KERNEL_MODE}")

add_library(
  rtekernels STATIC ${SOURCES} ${CMAKE_SOURCE_DIR}/rte-frontend/mo_rte_kind.F90
)
