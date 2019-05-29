lib_LIBRARIES = librrtmgp.a

librrtmgp_a_SOURCES =         \
    mo_gas_concentrations.F90 \
    mo_gas_optics.F90         \
    mo_gas_optics_rrtmgp.F90  \
    mo_rrtmgp_constants.F90   \
    mo_util_reorder.F90       \
    mo_util_string.F90

if ENABLE_OPENACC
kernel_dir = kernels-openacc
kernel_dir_extra = kernels
else
kernel_dir = kernels
kernel_dir_extra = kernels-openacc
endif

librrtmgp_a_SOURCES +=                      \
    $(kernel_dir)/mo_gas_optics_kernels.F90 \
    $(kernel_dir)/mo_reorder_kernels.F90

EXTRA_DIST =                                      \
    $(kernel_dir_extra)/mo_gas_optics_kernels.F90 \
    $(kernel_dir_extra)/mo_reorder_kernels.F90

dist_pkgdata_DATA =                        \
    data/rrtmgp-data-lw-g256-2018-12-04.nc \
    data/rrtmgp-data-sw-g224-2018-12-04.nc

if FCMODUC
mo_gas_concentrations_mod = MO_GAS_CONCENTRATIONS.$(FCMODEXT)
mo_gas_optics_kernels_mod = MO_GAS_OPTICS_KERNELS.$(FCMODEXT)
mo_gas_optics_mod = MO_GAS_OPTICS.$(FCMODEXT)
mo_gas_optics_rrtmgp_mod = MO_GAS_OPTICS_RRTMGP.$(FCMODEXT)
mo_reorder_kernels_mod = MO_REORDER_KERNELS.$(FCMODEXT)
mo_rrtmgp_constants_mod = MO_RRTMGP_CONSTANTS.$(FCMODEXT)
mo_util_reorder_mod = MO_UTIL_REORDER.$(FCMODEXT)
mo_util_string_mod = MO_UTIL_STRING.$(FCMODEXT)
else
mo_gas_concentrations_mod = mo_gas_concentrations.$(FCMODEXT)
mo_gas_optics_kernels_mod = mo_gas_optics_kernels.$(FCMODEXT)
mo_gas_optics_mod = mo_gas_optics.$(FCMODEXT)
mo_gas_optics_rrtmgp_mod = mo_gas_optics_rrtmgp.$(FCMODEXT)
mo_reorder_kernels_mod = mo_reorder_kernels.$(FCMODEXT)
mo_rrtmgp_constants_mod = mo_rrtmgp_constants.$(FCMODEXT)
mo_util_reorder_mod = mo_util_reorder.$(FCMODEXT)
mo_util_string_mod = mo_util_string.$(FCMODEXT)
endif

all_mods =                       \
    $(mo_gas_concentrations_mod) \
    $(mo_gas_optics_kernels_mod) \
    $(mo_gas_optics_mod)         \
    $(mo_gas_optics_rrtmgp_mod)  \
    $(mo_reorder_kernels_mod)    \
    $(mo_rrtmgp_constants_mod)   \
    $(mo_util_reorder_mod)       \
    $(mo_util_string_mod)

nodist_include_HEADERS = $(all_mods)

$(mo_gas_concentrations_mod): mo_gas_concentrations.$(OBJEXT)
mo_gas_concentrations.$(OBJEXT): $(mo_util_string_mod)

$(mo_gas_optics_kernels_mod): $(kernel_dir)/mo_gas_optics_kernels.$(OBJEXT)
$(kernel_dir)/mo_gas_optics_kernels.$(OBJEXT): $(mo_util_string_mod)

$(mo_gas_optics_mod): mo_gas_optics.$(OBJEXT)
mo_gas_optics.$(OBJEXT): $(mo_gas_concentrations_mod)

$(mo_gas_optics_rrtmgp_mod): mo_gas_optics_rrtmgp.$(OBJEXT)
mo_gas_optics_rrtmgp.$(OBJEXT):  \
    $(mo_gas_concentrations_mod) \
    $(mo_gas_optics_kernels_mod) \
    $(mo_gas_optics_mod)         \
    $(mo_rrtmgp_constants_mod)   \
    $(mo_util_reorder_mod)       \
    $(mo_util_string_mod)

$(mo_reorder_kernels_mod): $(kernel_dir)/mo_reorder_kernels.$(OBJEXT)

$(mo_rrtmgp_constants_mod): mo_rrtmgp_constants.$(OBJEXT)

$(mo_util_reorder_mod): mo_util_reorder.$(OBJEXT)
mo_util_reorder.$(OBJEXT): $(mo_reorder_kernels_mod)

$(mo_util_string_mod): mo_util_string.$(OBJEXT)

AM_FCFLAGS = $(FCMODOUT). $(FCMODINC)$(top_builddir)/rte

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

