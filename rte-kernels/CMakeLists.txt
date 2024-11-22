file(GLOB sourcefiles
    "mo_rte_solver_kernels.F90"
    "mo_fluxes_broadband_kernels.F90"
    "mo_rte_util_array.F90"
    "mo_optical_props_kernels.F90")

add_library(rtekernels STATIC ${sourcefiles})
