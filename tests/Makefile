#
# Location of RTE+RRTMGP libraries, module files.
#
RRTMGP_BUILD = $(RRTMGP_ROOT)/build
#
# RRTMGP library, module files
#
# LDFLAGS   += -L$(RRTMGP_BUILD)
# LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_BUILD)
#
# netcdf Fortran module files has to be in the search path or added via environment variable FCINCLUDE e.g. 
#FCINCLUDE += -I$(NFHOME)/include

# netcdf C and Fortran libraries have to be in the search path or added via environment variable LDFLAGS e.g. 
#LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LDFLAGS   += -L$(RRTMGP_BUILD)
LIBS      += -lrte -lrrtmgp
LIBS      += -lnetcdff -lnetcdf

VPATH  = .:$(RRTMGP_ROOT)/examples:$(RRTMGP_ROOT)/examples/rfmip-clear-sky:$(RRTMGP_ROOT)/examples/all-sky
VPATH += $(RRTMGP_ROOT)/rrtmgp-frontend:$(RRTMGP_ROOT)/extensions:$(RRTMGP_ROOT)/:$(RRTMGP_ROOT)/extensions/solar_variability

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)


#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_heating_rates.o mo_compute_bc.o mo_rrtmgp_clr_all_sky.o
ADDITIONS += mo_gas_optics_defs_rrtmgp.o
# File I/O
ADDITIONS += mo_simple_netcdf.o mo_rfmip_io.o
ADDITIONS += mo_testing_io.o mo_testing_utils.o 
# Cloud optics
CLOUDS    += mo_cloud_sampling.o mo_cloud_optics_rrtmgp.o mo_load_cloud_coefficients.o mo_garand_atmos_io.o
# Solar variability
ADDITIONS += mo_solar_variability.o

# Many codes will need to be updated if the library changes
# LIB_DEPS = $(RRTMGP_BUILD)/librte.a $(RRTMGP_BUILD)/librrtmgp.a
#
# Targets
#
all: check_variants check_equivalence test_zenith_angle_spherical_correction rte_sw_solver_unit_tests rte_optic_prop_unit_tests rte_lw_solver_unit_tests

check_equivalence:   $(ADDITIONS) $(LIB_DEPS) check_equivalence.o
check_equivalence.o: $(ADDITIONS) $(LIB_DEPS) check_equivalence.F90

check_variants:   $(ADDITIONS) $(LIB_DEPS) check_variants.o
check_variants.o: $(ADDITIONS) $(LIB_DEPS) check_variants.F90

test_zenith_angle_spherical_correction:   mo_zenith_angle_spherical_correction.o mo_rcemip_profiles.o $(ADDITIONS) $(LIB_DEPS) test_zenith_angle_spherical_correction.o
test_zenith_angle_spherical_correction.o: mo_zenith_angle_spherical_correction.o mo_rcemip_profiles.o $(ADDITIONS) $(LIB_DEPS) test_zenith_angle_spherical_correction.F90

mo_testing_io.o:        $(LIB_DEPS) mo_simple_netcdf.o mo_testing_io.F90

mo_cloud_optics_rrtmgp.o:     $(LIB_DEPS) mo_cloud_optics_rrtmgp.F90
mo_load_cloud_coefficients.o: $(LIB_DEPS) mo_simple_netcdf.o mo_cloud_optics_rrtmgp.o mo_load_cloud_coefficients.F90
mo_cloud_sampling.o:          $(LIB_DEPS) mo_cloud_sampling.F90

mo_gas_optics_defs_rrtmgp.o: $(LIB_DEPS) mo_testing_utils.o mo_simple_netcdf.o mo_gas_optics_defs_rrtmgp.F90

mo_load_coefficients.o: $(LIB_DEPS) mo_simple_netcdf.o mo_load_coefficients.F90
mo_rfmip_io.o.o:        $(LIB_DEPS) mo_simple_netcdf.o mo_rfmip_io.F90
mo_simple_netcdf.o:     $(LIB_DEPS) mo_simple_netcdf.F90

rte_optic_prop_unit_tests.o: $(LIB_DEPS) mo_testing_utils.o rte_optic_prop_unit_tests.F90
rte_optic_prop_unit_tests  : $(LIB_DEPS) mo_testing_utils.o rte_optic_prop_unit_tests.o

rte_lw_solver_unit_tests.o: $(LIB_DEPS) mo_testing_utils.o rte_lw_solver_unit_tests.F90
rte_lw_solver_unit_tests  : $(LIB_DEPS) mo_testing_utils.o rte_lw_solver_unit_tests.o

rte_sw_solver_unit_tests.o: $(LIB_DEPS) mo_testing_utils.o rte_sw_solver_unit_tests.F90
rte_sw_solver_unit_tests  : $(LIB_DEPS) mo_testing_utils.o rte_sw_solver_unit_tests.o


.PHONY: tests
tests:
	cp ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ./test_atmospheres.nc
	$(RUN_CMD) bash all_tests.sh


check:
	echo "Nothing to check in tests/"

clean:
	-rm clear_sky_regression *.o *.optrpt *.mod
