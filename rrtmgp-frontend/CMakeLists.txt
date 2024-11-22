file(GLOB sourcefiles
    "../gas-optics/mo_gas_optics_util_string.F90"
    "../gas-optics/mo_gas_concentrations.F90"
    "../gas-optics/mo_gas_optics.F90"
    "../gas-optics/mo_gas_optics_constants.F90"
    "mo_aerosol_optics_rrtmgp_merra.F90"
    "mo_gas_optics_rrtmgp.F90"
    "mo_cloud_optics_rrtmgp.F90")

add_library(rrtmgpf STATIC ${sourcefiles})
