message(STATUS "Configuring lib rrtmgpkernels")

# List of default files to compile. If KERNEL_MODE is accel or extern, then we
# try to compile accel versions or extern. If accel/extern versions of the files
# are not available, we fall back to the default
set(RTE_FORTRAN_KERNELS
    mo_cloud_optics_rrtmgp_kernels.F90 mo_gas_optics_rrtmgp_kernels.F90
)

find_files_in_folders(SOURCES "${RTE_FORTRAN_KERNELS}" "${KERNEL_MODE}")

add_library(rrtmgpkernels STATIC ${SOURCES})
target_link_libraries(rrtmgpkernels PRIVATE rtekernels)
