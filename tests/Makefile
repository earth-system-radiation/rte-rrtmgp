#
# Location of RTE+RRTMGP libraries, module files.
#   Default values can be overridden with environmental variables 
#
RRTMGP_DIR ?= $(RRTMGP_ROOT)/build
KERNEL_DIR ?= $(RRTMGP_ROOT)/build 

FCINCLUDE += -I$(RRTMGP_DIR)
LDFLAGS   += -L$(RRTMGP_DIR) -L$(KERNEL_DIR)
LIBS      += -lrrtmgpf -lrtef -lrtekernels -lrrtmgpkernels
#
# netcdf library, module files
# Environment variable NFHOME points to root of the Fortran interfaces 
#
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib
LIBS      += -lnetcdff
#
# Include libraries needed for CUDA kernels 
#
ifeq ($(RTE_KERNELS), extern)
	LIBS += -c++libs -cuda
endif

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
all: check_variants check_equivalence test_zenith_angle_spherical_correction

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

mo_load_coefficients.o: $(LIB_DEPS) mo_simple_netcdf.o mo_load_coefficients.F90
mo_rfmip_io.o.o:        $(LIB_DEPS) mo_simple_netcdf.o mo_rfmip_io.F90
mo_simple_netcdf.o:     $(LIB_DEPS) mo_simple_netcdf.F90

.PHONY: tests
tests:
	cp ${RRTMGP_DATA}/examples/rfmip-clear-sky/inputs/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ./test_atmospheres.nc
	$(RUN_CMD) bash all_tests.sh
check:
	echo "Nothing to check in tests/"

clean:
	-rm clear_sky_regression *.o *.optrpt *.mod
