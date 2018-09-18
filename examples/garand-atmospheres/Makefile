# Location of RTE+RRTMGP libraries, module files.
RRTMGP_ROOT = ../..
RRTMGP_BUILD = $(RRTMGP_ROOT)/build
# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
include $(RRTMGP_BUILD)/Makefile.conf
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_BUILD)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_BUILD)

#
# netcdf library, module files
# C and Fortran interfaces respectively
#
NCHOME = /opt/local
NFHOME = /Users/robert/Applications/$(FC)
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

VPATH = ./:$(RRTMGP_BUILD):$(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/examples

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)


#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_load_coefficients.o mo_simple_netcdf.o
ADDITIONS += mo_garand_atmos_io.o mo_heating_rates.o mo_fluxes_byband.o mo_fluxes_byband_kernels.o

#
# Targets
#
all: rrtmgp_garand_atmos

rrtmgp_garand_atmos: $(ADDITIONS) rrtmgp_garand_atmos.o

rrtmgp_garand_atmos.o: $(ADDITIONS) rrtmgp_garand_atmos.F90

mo_fluxes_byband.o: mo_fluxes_byband.F90 mo_fluxes_byband_kernels.o

mo_load_coefficients.o: mo_simple_netcdf.o mo_load_coefficients.F90
mo_garand_atmos_io.o:   mo_simple_netcdf.o mo_garand_atmos_io.F90

clean:
	-rm ../garand-atmos *.o *.optrpt ../*.optrpt *.mod
