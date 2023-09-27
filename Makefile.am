AM_FCFLAGS =
if USE_CBOOL
AM_FCFLAGS += $(FC_PP_DEF)RTE_USE_CBOOL
endif
if SP_ENABLED
AM_FCFLAGS += $(FC_PP_DEF)RTE_USE_SP
endif SP_ENABLED

LDADD = $(RPATH_FCFLAGS)

lib_LIBRARIES = librte.a

librte_a_SOURCES = \
  rte-frontend/mo_fluxes.F90 \
  rte-frontend/mo_optical_props.F90 \
  rte-frontend/mo_rte_config.F90 \
  rte-frontend/mo_rte_kind.F90 \
  rte-frontend/mo_rte_lw.F90 \
  rte-frontend/mo_rte_sw.F90 \
  rte-frontend/mo_rte_util_array_validation.F90 \
  rte-frontend/mo_source_functions.F90 \
  rte-kernels/mo_fluxes_broadband_kernels.F90 \
  rte-kernels/mo_rte_util_array.F90

if GPU_ENABLED
librte_a_SOURCES += \
  rte-kernels/accel/mo_optical_props_kernels.F90 \
  rte-kernels/accel/mo_rte_solver_kernels.F90
else !GPU_ENABLED
librte_a_SOURCES += \
  rte-kernels/mo_optical_props_kernels.F90 \
  rte-kernels/mo_rte_solver_kernels.F90
endif !GPU_ENABLED

lib_LIBRARIES += librrtmgp.a

librrtmgp_a_SOURCES = \
  gas-optics/mo_gas_concentrations.F90 \
  gas-optics/mo_gas_optics.F90 \
  gas-optics/mo_gas_optics_constants.F90 \
  gas-optics/mo_gas_optics_util_string.F90 \
  rrtmgp-frontend/mo_gas_optics_rrtmgp.F90

if GPU_ENABLED
librrtmgp_a_SOURCES += \
  rrtmgp-kernels/accel/mo_gas_optics_rrtmgp_kernels.F90
else !GPU_ENABLED
librrtmgp_a_SOURCES += \
  rrtmgp-kernels/mo_gas_optics_rrtmgp_kernels.F90
endif !GPU_ENABLED

if TESTS_ENABLED

AM_FCFLAGS += $(NETCDF_FCFLAGS)
LDADD += librrtmgp.a librte.a $(NETCDF_FCLIBS)

check_PROGRAMS = \
  examples/all-sky/rrtmgp_allsky \
  examples/rfmip-clear-sky/rrtmgp_rfmip_lw \
  examples/rfmip-clear-sky/rrtmgp_rfmip_sw \
  tests/check_equivalence \
  tests/check_variants \
  tests/test_zenith_angle_spherical_correction

examples_all_sky_rrtmgp_allsky_SOURCES = \
  examples/all-sky/mo_load_aerosol_coefficients.F90 \
  examples/all-sky/mo_load_cloud_coefficients.F90 \
  examples/all-sky/rrtmgp_allsky.F90 \
  examples/mo_load_coefficients.F90 \
  examples/mo_simple_netcdf.F90 \
  rrtmgp-frontend/mo_aerosol_optics_rrtmgp_merra.F90 \
  rrtmgp-frontend/mo_cloud_optics_rrtmgp.F90

examples_rfmip_clear_sky_rrtmgp_rfmip_lw_SOURCES = \
  examples/mo_load_coefficients.F90 \
  examples/mo_simple_netcdf.F90 \
  examples/rfmip-clear-sky/mo_rfmip_io.F90 \
  examples/rfmip-clear-sky/rrtmgp_rfmip_lw.F90

examples_rfmip_clear_sky_rrtmgp_rfmip_sw_SOURCES = \
  examples/mo_load_coefficients.F90 \
  examples/mo_simple_netcdf.F90 \
  examples/rfmip-clear-sky/mo_rfmip_io.F90 \
  examples/rfmip-clear-sky/rrtmgp_rfmip_sw.F90

tests_check_equivalence_SOURCES = \
  examples/mo_load_coefficients.F90 \
  examples/mo_simple_netcdf.F90 \
  examples/rfmip-clear-sky/mo_rfmip_io.F90 \
  extensions/mo_heating_rates.F90 \
  tests/check_equivalence.F90

tests_check_variants_SOURCES = \
  examples/mo_load_coefficients.F90 \
  examples/mo_simple_netcdf.F90 \
  examples/rfmip-clear-sky/mo_rfmip_io.F90 \
  extensions/mo_heating_rates.F90 \
  tests/check_variants.F90 \
  tests/mo_testing_io.F90

tests_test_zenith_angle_spherical_correction_SOURCES = \
  examples/mo_load_coefficients.F90 \
  examples/mo_simple_netcdf.F90 \
  extensions/mo_heating_rates.F90 \
  extensions/mo_zenith_angle_spherical_correction.F90 \
  tests/mo_rcemip_profiles.F90 \
  tests/test_zenith_angle_spherical_correction.F90

endif TESTS_ENABLED

TESTS = \
  examples/all-sky/all-sky.test \
  examples/rfmip-clear-sky/rfmip-clear-sky.test \
  tests/check_equivalence.test \
  tests/check_variants.test \
  tests/test_zenith_angle_spherical_correction.test

$(TEST_LOGS): @TESTS_ENABLED_TRUE@ $(testcachedir)/$(am__dirstamp)
TEST_DATA_URL= https://github.com/earth-system-radiation/rrtmgp-data/archive
TEST_DATA_TAG= df02975ab93165b34a59f0d04b4ae6148fe5127c
$(testcachedir)/$(am__dirstamp):
	@$(MKDIR_P) $(testcachedir) && $(am__cd) $(testcachedir) && \
	echo "Fetching test input data..." && \
	$(CURL) -sL $(TEST_DATA_URL)/$(TEST_DATA_TAG).tar.gz | \
	  $(AMTAR) xz --strip-components=1 && touch $(am__dirstamp)

mostlyclean-local:
	-rm -f *.$(FC_MOD_FILE_EXT)

clean-local:
	-rm -f examples/all-sky/*.nc examples/rfmip-clear-sky/*.nc tests/*.nc

distclean-local:
	-rm -rf Makefile.exe.deps Makefile.lib.deps $(testcachedir)
	-rm -rf mkhelper/depgen/*.pyc mkhelper/depgen/*.pyo mkhelper/depgen/__pycache__

# Compile *.F90 without CPPFLAGS, which are normally meant for the C compiler
# and might not be compatible with the Fortran compiler:
PPFCCOMPILE = $(FC) $(AM_FCFLAGS) $(FCFLAGS)

# Handle Fortran dependency files:
./Makefile.exe.deps: Makefile \
  $(examples_all_sky_rrtmgp_allsky_SOURCES) \
  $(examples_rfmip_clear_sky_rrtmgp_rfmip_lw_SOURCES) \
  $(examples_rfmip_clear_sky_rrtmgp_rfmip_sw_SOURCES) \
  $(tests_check_equivalence_SOURCES) \
  $(tests_check_variants_SOURCES) \
  $(tests_test_zenith_angle_spherical_correction_SOURCES)

./Makefile.lib.deps: Makefile \
  $(librrtmgp_a_SOURCES) \
  $(librte_a_SOURCES)

./Makefile.exe.deps ./Makefile.lib.deps:
	$(AM_V_GEN)list=`echo $^ | tr ' ' '\n' | grep '\.F90$$' | $(am__uniquify_input)`; \
	$(am__nobase_strip_setup); obj_names=`$(am__nobase_strip) | sed 's/.F90$$/.@OBJEXT@/g'`; \
	$(PYTHON) $(top_srcdir)/mkhelper/depgen.py --src-root='$(top_srcdir)' \
	  --pp-enable --pp-eval-expr --pp-macro-flag='$(FC_PP_DEF)' \
	  --pp-inc-flag='$(FC_INC_FLAG_PP_f90)' \
	  --pp-inc-order='$(FC_INC_ORDER_PP_f90)' --fc-enable \
	  --fc-mod-ext='$(FC_MOD_FILE_EXT)' \
	  --fc-mod-upper='$(FC_MOD_FILE_UPPER)' --fc-inc-flag='$(FC_INC_FLAG)' \
	  --fc-inc-order='$(FC_INC_ORDER)' --fc-mod-dir-flag='$(FC_MOD_OUT)' \
	  --fc-external-mods='netcdf' --obj-name $$obj_names --input $$list \
	  -- $(DEPGEN_FCFLAGS) $(AM_FCFLAGS) $(FCFLAGS) > $@~ && mv $@~ $@

install-data-local: $(lib_LIBRARIES)
	@$(NORMAL_INSTALL)
	@list=`$(PYTHON) $(top_srcdir)/mkhelper/deplist.py -f ./Makefile.lib.deps | sed -n '/.$(FC_MOD_FILE_EXT)$$/p'`; \
	test -n "$(includedir)" || list=; \
	if test -n "$$list"; then \
	  echo " $(MKDIR_P) '$(DESTDIR)$(includedir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(includedir)" || exit 1; \
	fi; \
	for p in $$list; do \
	  echo " $(INSTALL_HEADER) $$p '$(DESTDIR)$(includedir)'"; \
	  $(INSTALL_HEADER) $$p "$(DESTDIR)$(includedir)" || exit $$?; \
	done

uninstall-local:
	@$(NORMAL_UNINSTALL)
	@list=`$(PYTHON) $(top_srcdir)/mkhelper/deplist.py -f ./Makefile.lib.deps | sed -n '/.$(FC_MOD_FILE_EXT)$$/p'`; \
	test -n "$(includedir)" || list=; \
	files=$$list; dir='$(DESTDIR)$(includedir)'; $(am__uninstall_files_from_dir)

EXTRA_DIST = \
  autogen.sh \
  examples/compare-to-reference.py \
  mkhelper/_argparse.py \
  mkhelper/depgen.py \
  mkhelper/depgen/__init__.py \
  mkhelper/depgen/fortran_parser.py \
  mkhelper/depgen/line_control.py \
  mkhelper/depgen/preprocessor.py \
  mkhelper/deplist.py \
  mkhelper/fortmodcmp.py

@TESTS_ENABLED_TRUE@@am__include@ @am__quote@./Makefile.exe.deps@am__quote@ # am--include-marker
@am__include@ @am__quote@./Makefile.lib.deps@am__quote@ # am--include-marker
