cmake_minimum_required(VERSION 3.0)

set(FORTRAN_FILES ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/mo_optical_props.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/mo_rte_sw.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/mo_rte_util_array.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/mo_source_functions.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/kernels/mo_optical_props_kernels.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/kernels/mo_fluxes_broadband_kernels.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/kernels/mo_rte_solver_kernels.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/mo_fluxes.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/mo_rte_lw.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rte/mo_rte_kind.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/mo_rrtmgp_constants.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/mo_gas_optics.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/mo_rrtmgp_util_string.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/mo_gas_concentrations.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/mo_gas_optics_rrtmgp.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/mo_rrtmgp_util_reorder.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../rrtmgp/kernels/mo_gas_optics_kernels.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../examples/mo_simple_netcdf.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../examples/all-sky/rrtmgp_allsky.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../examples/all-sky/mo_garand_atmos_io.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../examples/all-sky/mo_load_cloud_coefficients.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../examples/mo_load_coefficients.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/mo_fluxes_byband.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/mo_fluxes_byband_kernels.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/cloud_optics/mo_cloud_sampling.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/cloud_optics/mo_cloud_optics.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/mo_compute_bc.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/mo_rrtmgp_clr_all_sky.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/mo_fluxes_bygpoint.F90
                  ${CMAKE_CURRENT_SOURCE_DIR}/../../../extensions/mo_heating_rates.F90
               )
add_executable(allsky_fortran ${FORTRAN_FILES})
set_source_files_properties(${FORTRAN_FILES} PROPERTIES COMPILE_FLAGS "${YAKL_F90_FLAGS}")
target_link_libraries(allsky_fortran "${F90_LINK}")

