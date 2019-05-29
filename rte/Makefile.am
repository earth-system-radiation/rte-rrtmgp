lib_LIBRARIES = librte.a

librte_a_SOURCES =          \
    mo_fluxes.F90           \
    mo_optical_props.F90    \
    mo_rte_kind.F90         \
    mo_rte_lw.F90           \
    mo_rte_sw.F90           \
    mo_source_functions.F90 \
    mo_util_array.F90

if ENABLE_OPENACC
kernel_dir = kernels-openacc
kernel_dir_extra = kernels
else
kernel_dir = kernels
kernel_dir_extra = kernels-openacc
endif

librte_a_SOURCES +=                               \
    $(kernel_dir)/mo_fluxes_broadband_kernels.F90 \
    $(kernel_dir)/mo_optical_props_kernels.F90    \
    $(kernel_dir)/mo_rte_solver_kernels.F90

EXTRA_DIST =                                            \
    $(kernel_dir_extra)/mo_fluxes_broadband_kernels.F90 \
    $(kernel_dir_extra)/mo_optical_props_kernels.F90    \
    $(kernel_dir_extra)/mo_rte_solver_kernels.F90

if FCMODUC
mo_fluxes_broadband_kernels_mod = MO_FLUXES_BROADBAND_KERNELS.$(FCMODEXT)
mo_fluxes_mod = MO_FLUXES.$(FCMODEXT)
mo_optical_props_kernels_mod = MO_OPTICAL_PROPS_KERNELS.$(FCMODEXT)
mo_optical_props_mod = MO_OPTICAL_PROPS.$(FCMODEXT)
mo_rte_kind_mod = MO_RTE_KIND.$(FCMODEXT)
mo_rte_lw_mod = MO_RTE_LW.$(FCMODEXT)
mo_rte_solver_kernels_mod = MO_RTE_SOLVER_KERNELS.$(FCMODEXT)
mo_rte_sw_mod = MO_RTE_SW.$(FCMODEXT)
mo_source_functions_mod = MO_SOURCE_FUNCTIONS.$(FCMODEXT)
mo_util_array_mod = MO_UTIL_ARRAY.$(FCMODEXT)
else
mo_fluxes_broadband_kernels_mod = mo_fluxes_broadband_kernels.$(FCMODEXT)
mo_fluxes_mod = mo_fluxes.$(FCMODEXT)
mo_optical_props_kernels_mod = mo_optical_props_kernels.$(FCMODEXT)
mo_optical_props_mod = mo_optical_props.$(FCMODEXT)
mo_rte_kind_mod = mo_rte_kind.$(FCMODEXT)
mo_rte_lw_mod = mo_rte_lw.$(FCMODEXT)
mo_rte_solver_kernels_mod = mo_rte_solver_kernels.$(FCMODEXT)
mo_rte_sw_mod = mo_rte_sw.$(FCMODEXT)
mo_source_functions_mod = mo_source_functions.$(FCMODEXT)
mo_util_array_mod = mo_util_array.$(FCMODEXT)
endif

all_mods =                             \
    $(mo_fluxes_broadband_kernels_mod) \
    $(mo_fluxes_mod)                   \
    $(mo_optical_props_kernels_mod)    \
    $(mo_optical_props_mod)            \
    $(mo_rte_kind_mod)                 \
    $(mo_rte_lw_mod)                   \
    $(mo_rte_solver_kernels_mod)       \
    $(mo_rte_sw_mod)                   \
    $(mo_source_functions_mod)         \
    $(mo_util_array_mod)

nodist_include_HEADERS = $(all_mods)

$(mo_fluxes_broadband_kernels_mod): $(kernel_dir)/mo_fluxes_broadband_kernels.$(OBJEXT)
$(kernel_dir)/mo_fluxes_broadband_kernels.$(OBJEXT): $(mo_rte_kind_mod)

$(mo_fluxes_mod): mo_fluxes.$(OBJEXT)
mo_fluxes.$(OBJEXT):                   \
    $(mo_fluxes_broadband_kernels_mod) \
    $(mo_optical_props_mod)            \
    $(mo_rte_kind_mod)

$(mo_optical_props_kernels_mod): $(kernel_dir)/mo_optical_props_kernels.$(OBJEXT)
$(kernel_dir)/mo_optical_props_kernels.$(OBJEXT): $(mo_rte_kind_mod)

$(mo_optical_props_mod): mo_optical_props.$(OBJEXT)
mo_optical_props.$(OBJEXT):         \
    $(mo_optical_props_kernels_mod) \
    $(mo_rte_kind_mod)              \
    $(mo_util_array_mod)

$(mo_rte_kind_mod): mo_rte_kind.$(OBJEXT)

$(mo_rte_lw_mod): mo_rte_lw.$(OBJEXT)
mo_rte_lw.$(OBJEXT):             \
    $(mo_fluxes_mod)             \
    $(mo_optical_props_mod)      \
    $(mo_rte_kind_mod)           \
    $(mo_rte_solver_kernels_mod) \
    $(mo_source_functions_mod)

$(mo_rte_solver_kernels_mod): $(kernel_dir)/mo_rte_solver_kernels.$(OBJEXT)
$(kernel_dir)/mo_rte_solver_kernels.$(OBJEXT): $(mo_rte_kind_mod)

$(mo_rte_sw_mod): mo_rte_sw.$(OBJEXT)
mo_rte_sw.$(OBJEXT):             \
    $(mo_fluxes_mod)             \
    $(mo_optical_props_mod)      \
    $(mo_rte_kind_mod)           \
    $(mo_rte_solver_kernels_mod)

$(mo_source_functions_mod): mo_source_functions.$(OBJEXT)
mo_source_functions.$(OBJEXT): \
    $(mo_optical_props_mod)    \
    $(mo_rte_kind_mod)

$(mo_util_array_mod): mo_util_array.$(OBJEXT)
mo_util_array.$(OBJEXT): $(mo_rte_kind_mod)

AM_FCFLAGS = $(FCMODOUT).

# Override the automatically generated rule
.F90.o:
	$(AM_V_PPFC)$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ $(FCFLAGS_F90) $<

# Regenerate a module file if it does not exist
$(all_mods):
	@if test -z $<; then \
          echo "Cannot find Fortran source file providing module '$(basename $(@F))'." >&2; \
        elif test ! -f $@; then \
          rm -f $<; $(MAKE) $(AM_MAKEFLAGS) $<; \
        fi

mostlyclean-local:
	rm -f *.i

MOSTLYCLEANFILES = $(all_mods)

