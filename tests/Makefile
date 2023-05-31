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
# netcdf library, module files
# Environment variables NCHOME and NFHOME point to root of C and Fortran interfaces respectively -
#
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(RRTMGP_BUILD)
LIBS      += -lrte -lrrtmgp
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
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
# File I/O
ADDITIONS += mo_load_coefficients.o mo_simple_netcdf.o mo_rfmip_io.o
ADDITIONS += mo_testing_io.o
# Cloud optics
CLOUDS    += mo_cloud_sampling.o mo_cloud_optics_rrtmgp.o mo_load_cloud_coefficients.o mo_garand_atmos_io.o
# Solar variability
ADDITIONS += mo_solar_variability.o

# Many codes will need to be updated if the library changes
# LIB_DEPS = $(RRTMGP_BUILD)/librte.a $(RRTMGP_BUILD)/librrtmgp.a
#
# Targets
#
all: clear_sky_regression test_zenith_angle_spherical_correction

clear_sky_regression:   $(ADDITIONS) $(LIB_DEPS) clear_sky_regression.o
clear_sky_regression.o: $(ADDITIONS) $(LIB_DEPS) clear_sky_regression.F90

test_zenith_angle_spherical_correction:   mo_zenith_angle_spherical_correction.o mo_rcemip_profiles.o $(ADDITIONS) $(LIB_DEPS) test_zenith_angle_spherical_correction.o
test_zenith_angle_spherical_correction.o: mo_zenith_angle_spherical_correction.o mo_rcemip_profiles.o $(ADDITIONS) $(LIB_DEPS) test_zenith_angle_spherical_correction.F90

mo_testing_io.o:        $(LIB_DEPS) mo_simple_netcdf.o mo_testing_io.F90

mo_cloud_optics_rrtmgp.o:     $(LIB_DEPS) mo_cloud_optics_rrtmgp.F90
mo_load_cloud_coefficients.o: $(LIB_DEPS) mo_simple_netcdf.o mo_cloud_optics_rrtmgp.o mo_load_cloud_coefficients.F90
mo_cloud_sampling.o:          $(LIB_DEPS) mo_cloud_sampling.F90

mo_load_coefficients.o: $(LIB_DEPS) mo_simple_netcdf.o mo_load_coefficients.F90
mo_rfmip_io.o.o:        $(LIB_DEPS) mo_simple_netcdf.o mo_rfmip_io.F90
mo_simple_netcdf.o:     $(LIB_DEPS) mo_simple_netcdf.F90

.PHONY: tests
tests:
	cp ${RRTMGP_ROOT}/examples/rfmip-clear-sky/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ./test_atmospheres.nc
	$(RUN_CMD) ./clear_sky_regression test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g128.nc
	$(RUN_CMD) ./clear_sky_regression test_atmospheres.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g112.nc

check:
	python verification.py

clean:
	-rm clear_sky_regression *.o *.optrpt *.mod
