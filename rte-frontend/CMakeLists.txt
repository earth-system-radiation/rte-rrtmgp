file(GLOB sourcefiles
    "mo_optical_props.F90"
    "mo_rte_lw.F90"
    "mo_fluxes.F90"
    "mo_rte_sw.F90"
    "mo_rte_util_array_validation.F90"
    "mo_rte_config.F90"
    "mo_source_functions.F90")

add_library(rtef STATIC ${sourcefiles})
