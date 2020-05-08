kernels = kernels

lib_LTLIBRARIES = librte.la librrtmgp.la

#RTE solver library.
rte = $(top_srcdir)/rte
librte_la_SOURCES = $(rte)/mo_fluxes.F90 \
                    $(rte)/mo_optical_props.F90 \
                    $(rte)/mo_rte_kind.F90 \
                    $(rte)/mo_rte_lw.F90 \
                    $(rte)/mo_rte_sw.F90 \
                    $(rte)/mo_rte_util_array.F90 \
                    $(rte)/mo_source_functions.F90 \
                    $(rte)/$(kernels)/mo_fluxes_broadband_kernels.F90 \
                    $(rte)/$(kernels)/mo_optical_props_kernels.F90 \
                    $(rte)/$(kernels)/mo_rte_solver_kernels.F90
$(rte)/mo_fluxes.lo: $(rte)/mo_rte_kind.lo \
                     $(rte)/$(kernels)/mo_fluxes_broadband_kernels.lo \
                     $(rte)/mo_optical_props.lo \
                     $(rte)/mo_rte_util_array.lo
$(rte)/mo_optical_props.lo: $(rte)/mo_rte_kind.lo \
                            $(rte)/mo_rte_util_array.lo \
                            $(rte)/$(kernels)/mo_optical_props_kernels.lo
$(rte)/mo_rte_lw.lo: $(rte)/mo_rte_kind.lo \
                     $(rte)/mo_rte_util_array.lo \
                     $(rte)/mo_optical_props.lo \
                     $(rte)/mo_source_functions.lo \
                     $(rte)/mo_fluxes.lo \
                     $(rte)/$(kernels)/mo_rte_solver_kernels.lo
$(rte)/mo_rte_sw.lo: $(rte)/mo_rte_kind.lo \
                     $(rte)/mo_rte_util_array.lo \
                     $(rte)/mo_optical_props.lo \
                     $(rte)/mo_source_functions.lo \
                     $(rte)/mo_fluxes.lo \
                     $(rte)/$(kernels)/mo_rte_solver_kernels.lo
$(rte)/mo_rte_util_array.lo: $(rte)/mo_rte_kind.lo
$(rte)/mo_source_functions.lo: $(rte)/mo_rte_kind.lo \
                               $(rte)/mo_optical_props.lo
$(rte)/$(kernels)/mo_fluxes_broadband_kernels.lo: $(rte)/mo_rte_kind.lo
$(rte)/$(kernels)/mo_optical_props_kernels.lo: $(rte)/mo_rte_kind.lo
$(rte)/$(kernels)/mo_rte_solver_kernels.lo: $(rte)/mo_rte_kind.lo


#RRTMGP gas optics library.
rrtmgp = $(top_srcdir)/rrtmgp
librrtmgp_la_SOURCES = $(rrtmgp)/mo_gas_concentrations.F90 \
                       $(rrtmgp)/mo_gas_optics.F90 \
                       $(rrtmgp)/mo_gas_optics_rrtmgp.F90 \
                       $(rrtmgp)/mo_rrtmgp_constants.F90 \
                       $(rrtmgp)/mo_rrtmgp_util_reorder.F90 \
                       $(rrtmgp)/mo_rrtmgp_util_string.F90 \
                       $(rrtmgp)/$(kernels)/mo_gas_optics_kernels.F90 \
                       $(rrtmgp)/$(kernels)/mo_rrtmgp_util_reorder_kernels.F90 \
                       $(top_srcdir)/examples/mo_load_coefficients.F90
librrtmgp_la_LIBADD = $(top_builddir)/librte.la
$(rrtmgp)/mo_gas_concentrations.lo: $(rrtmgp)/mo_rrtmgp_util_string.lo
$(rrtmgp)/mo_gas_optics.lo: $(rrtmgp)/mo_gas_concentrations.lo
$(rrtmgp)/mo_gas_optics_rrtmgp.lo: $(rrtmgp)/mo_rrtmgp_constants.lo \
                                   $(rrtmgp)/mo_rrtmgp_util_reorder.lo \
                                   $(rrtmgp)/mo_rrtmgp_util_string.lo \
                                   $(rrtmgp)/$(kernels)/mo_gas_optics_kernels.lo \
                                   $(rrtmgp)/mo_gas_concentrations.lo \
                                   $(rrtmgp)/mo_gas_optics.lo
$(rrtmgp)/mo_rrtmgp_util_reorder.lo: $(rrtmgp)/$(kernels)/mo_rrtmgp_util_reorder_kernels.lo

include_HEADERS = mo_fluxes.mod \
                  mo_optical_props.mod \
                  mo_rte_kind.mod \
                  mo_rte_lw.mod \
                  mo_rte_sw.mod \
                  mo_rte_util_array.mod \
                  mo_source_functions.mod \
                  mo_fluxes_broadband_kernels.mod \
                  mo_optical_props_kernels.mod \
                  mo_rte_solver_kernels.mod \
                  mo_gas_concentrations.mod \
                  mo_gas_optics.mod \
                  mo_gas_optics_rrtmgp.mod \
                  mo_rrtmgp_constants.mod \
                  mo_rrtmgp_util_reorder.mod \
                  mo_rrtmgp_util_string.mod \
                  mo_gas_optics_kernels.mod \
                  mo_rrtmgp_util_reorder_kernels.mod \
                  mo_load_coefficients.mod

CLEANFILES = *.mod
