#
# Location of RTE+RRTMGP libraries, module files.
#
RRTMGP_ROOT  = ../..
NFHOME       = /usr
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

VPATH = ../:$(RRTMGP_ROOT)/extensions/cloud_optics

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_load_coefficients.o mo_simple_netcdf.o mo_cloud_optics.o mo_load_cloud_coefficients.o
ADDITIONS += mo_garand_atmos_io.o

#
# Targets
#
all: rrtmgp_allsky

rrtmgp_allsky: $(ADDITIONS) rrtmgp_allsky.o

rrtmgp_allsky.o: $(ADDITIONS) rrtmgp_allsky.F90

mo_cloud_optics.o: mo_cloud_optics.F90
mo_load_coefficients.o: mo_simple_netcdf.o mo_load_coefficients.F90
mo_garand_atmos_io.o:   mo_simple_netcdf.o mo_garand_atmos_io.F90
mo_load_cloud_coefficients.o: mo_simple_netcdf.o mo_cloud_optics.o mo_load_cloud_coefficients.F90

tests:
	cp garand-atmos-1.nc rrtmgp-allsky.nc
	$(RUN_CMD) ./rrtmgp_allsky rrtmgp-allsky.nc ${RRTMGP_ROOT}/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc ${RRTMGP_ROOT}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc 128
	$(RUN_CMD) ./rrtmgp_allsky rrtmgp-allsky.nc ${RRTMGP_ROOT}/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc ${RRTMGP_ROOT}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc 128

check:
	python3 ./compare-to-reference.py

clean:
	-rm rrtmgp_allsky *.o *.optrpt ../*.optrpt *.mod
