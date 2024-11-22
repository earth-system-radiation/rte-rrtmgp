file(GLOB sourcefiles
    # "../rte-frontend/mo_rte_kind.F90"
    "mo_gas_optics_rrtmgp_kernels.F90")

add_library(rrtmgpkernels STATIC ${sourcefiles})
