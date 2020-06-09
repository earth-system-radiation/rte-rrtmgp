#
# Location of RTE+RRTMGP libraries, module files.
#
RRTMGP_BUILD = $(RRTMGP_ROOT)/build
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_BUILD)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_BUILD)
#
# netcdf library, module files
# Environment variables NCHOME and NFHOME point to root of C and Fortran interfaces respectively -
#
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

VPATH  = .:$(RRTMGP_ROOT)/examples:$(RRTMGP_ROOT)/examples/rfmip-clear-sky
VPATH += $(RRTMGP_ROOT)/extensions:$(RRTMGP_ROOT)/extensions/cloud_optics:$(RRTMGP_ROOT)/extensions/solar_variability

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
CLOUDS    += mo_cloud_sampling.o mo_cloud_optics.o mo_load_cloud_coefficients.o
# Solar variability
ADDITIONS += mo_solar_variability.o

# Many codes will need to be updated if the library changes
LIB_DEPS = $(RRTMGP_BUILD)/librte.a $(RRTMGP_BUILD)/librrtmgp.a
#
# Targets
#
all: clear_sky_regression

clear_sky_regression:   $(ADDITIONS) $(LIB_DEPS) clear_sky_regression.o
clear_sky_regression.o: $(ADDITIONS) $(LIB_DEPS) clear_sky_regression.F90

mo_testing_io.o:        $(LIB_DEPS) mo_simple_netcdf.o mo_testing_io.F90

mo_cloud_optics.o:            $(LIB_DEPS) mo_cloud_optics.F90
mo_load_cloud_coefficients.o: $(LIB_DEPS) mo_simple_netcdf.o mo_cloud_optics.o mo_load_cloud_coefficients.F90
mo_cloud_sampling.o:          $(LIB_DEPS) mo_cloud_sampling.F90

mo_load_coefficients.o: $(LIB_DEPS) mo_simple_netcdf.o mo_load_coefficients.F90
mo_rfmip_io.o.o:        $(LIB_DEPS) mo_simple_netcdf.o mo_rfmip_io.F90
mo_simple_netcdf.o:     $(LIB_DEPS) mo_simple_netcdf.F90

clean:
	-rm clear_sky_regression *.o *.optrpt
