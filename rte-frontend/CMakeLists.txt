file(GLOB sourcefiles
    "mo_rte_kind.F90"
    "mo_optical_props.F90"
    "mo_rte_lw.F90"
    "mo_fluxes.F90"
    "mo_rte_sw.F90"
    "mo_rte_util_array_validation.F90"
    "mo_rte_config.F90"
    "mo_source_functions.F90")

add_library(rte-frontend STATIC ${sourcefiles})
# target_include_directories(microhh_core PRIVATE "../include" "../include_kernels" "../extern/pybind11/include" SYSTEM ${INCLUDE_DIRS})
# target_link_libraries(microhh_core PRIVATE pybind11::embed ${LIBS})
