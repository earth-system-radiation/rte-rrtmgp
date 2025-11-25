add_library(ssm OBJECT)

target_sources(
    ssm
    PRIVATE # cmake-format: sort
    mo_gas_optics_ssm.F90
)