check_PROGRAMS =

if HAVE_NETCDF
check_PROGRAMS +=                   \
    rfmip-clear-sky/rrtmgp_rfmip_lw \
    rfmip-clear-sky/rrtmgp_rfmip_sw
endif

TESTS = rfmip-clear-sky.test

rfmip_clear_sky_rrtmgp_rfmip_lw_SOURCES = \
    mo_load_coefficients.F90              \
    mo_simple_netcdf.F90                  \
    rfmip-clear-sky/mo_rfmip_io.F90       \
    rfmip-clear-sky/rrtmgp_rfmip_lw.F90

rfmip_clear_sky_rrtmgp_rfmip_sw_SOURCES = \
    mo_load_coefficients.F90              \
    mo_simple_netcdf.F90                  \
    rfmip-clear-sky/mo_rfmip_io.F90       \
    rfmip-clear-sky/rrtmgp_rfmip_sw.F90

EXTRA_DIST =                                \
    rfmip-clear-sky/compare-to-reference.py \
    rfmip-clear-sky/run-rfmip-examples.py   \
    rfmip-clear-sky/stage_files.py

if FCMODUC
mo_load_coefficients_mod = MO_LOAD_COEFFICIENTS.$(FCMODEXT)
mo_rfmip_io_mod = MO_RFMIP_IO.$(FCMODEXT)
mo_simple_netcdf_mod = MO_SIMPLE_NETCDF.$(FCMODEXT)
else
mo_load_coefficients_mod = mo_load_coefficients.$(FCMODEXT)
mo_rfmip_io_mod = mo_rfmip_io.$(FCMODEXT)
mo_simple_netcdf_mod = mo_simple_netcdf.$(FCMODEXT)
endif

all_mods =                      \
    $(mo_load_coefficients_mod) \
    $(mo_rfmip_io_mod)          \
    $(mo_simple_netcdf_mod)

$(mo_load_coefficients_mod): mo_load_coefficients.$(OBJEXT)
mo_load_coefficients.$(OBJEXT): $(mo_simple_netcdf_mod)

$(mo_rfmip_io_mod): rfmip-clear-sky/mo_rfmip_io.$(OBJEXT)
rfmip-clear-sky/mo_rfmip_io.$(OBJEXT): $(mo_simple_netcdf_mod)

$(mo_simple_netcdf_mod): mo_simple_netcdf.$(OBJEXT)

rfmip-clear-sky/rrtmgp_rfmip_lw.$(OBJEXT): $(mo_load_coefficients_mod) $(mo_rfmip_io_mod)

rfmip-clear-sky/rrtmgp_rfmip_sw.$(OBJEXT): $(mo_load_coefficients_mod) $(mo_rfmip_io_mod)

AM_FCFLAGS = $(FCMODOUT). $(FCMODINC)$(top_builddir)/rte $(FCMODINC)$(top_builddir)/rrtmgp $(NETCDF_FCFLAGS) $(GPTL_FCFLAGS)
LDADD = $(top_builddir)/rrtmgp/librrtmgp.a $(top_builddir)/rte/librte.a $(GPTL_FCLIBS) $(NETCDF_FCLIBS)

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

CLEANFILES = \
    rfmip-clear-sky/generate-output-file-templates.py                               \
    rfmip-clear-sky/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc \
    rfmip-clear-sky/reference/rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc      \
    rfmip-clear-sky/reference/rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc      \
    rfmip-clear-sky/reference/rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc      \
    rfmip-clear-sky/reference/rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc      \
    rfmip-clear-sky/rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc                \
    rfmip-clear-sky/rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc                \
    rfmip-clear-sky/rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc                \
    rfmip-clear-sky/rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc                \
    rfmip-clear-sky/timing.8

AUTOMAKE_OPTIONS = color-tests serial-tests

