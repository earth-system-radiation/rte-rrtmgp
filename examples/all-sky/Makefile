#
# Location of RTE+RRTMGP libraries, module files.
#   Default values can be overridden with environmental variables 
#
RRTMGP_DIR ?= $(RRTMGP_ROOT)/build
KERNEL_DIR ?= $(RRTMGP_ROOT)/build 
#
# RRTMGP library, module files
#
FCINCLUDE += -I$(RRTMGP_DIR)
LDFLAGS   += -L$(RRTMGP_DIR) -L$(KERNEL_DIR)
LIBS      += -lrtef -lrrtmgpf -lrtekernels -lrrtmgpkernels

#
# netcdf library, module files
# Environment variable NFHOME points to root of the Fortran interface
#
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib
LIBS      += -lnetcdff

VPATH = ../:$(RRTMGP_ROOT)/rrtmgp-frontend # Needed for cloud_optics and aerosol_optics

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_load_coefficients.o mo_simple_netcdf.o mo_cloud_optics_rrtmgp.o mo_load_cloud_coefficients.o
ADDITIONS += mo_garand_atmos_io.o
ADDITIONS += mo_aerosol_optics_rrtmgp_merra.o mo_load_aerosol_coefficients.o

#
# Targets
#
all: rrtmgp_allsky

rrtmgp_allsky: $(ADDITIONS) rrtmgp_allsky.o

rrtmgp_allsky.o: $(ADDITIONS) rrtmgp_allsky.F90

mo_cloud_optics_rrtmgp.o: mo_cloud_optics_rrtmgp.F90
mo_load_coefficients.o:   mo_simple_netcdf.o mo_load_coefficients.F90
mo_garand_atmos_io.o:     mo_simple_netcdf.o mo_garand_atmos_io.F90
mo_load_cloud_coefficients.o: mo_simple_netcdf.o mo_cloud_optics_rrtmgp.o mo_load_cloud_coefficients.F90
mo_aerosol_optics_rrtmgp_merra.o: mo_aerosol_optics_rrtmgp_merra.F90
mo_load_aerosol_coefficients.o: mo_simple_netcdf.o mo_aerosol_optics_rrtmgp_merra.o mo_load_aerosol_coefficients.F90

tests:
	cp ${RRTMGP_DATA}/examples/all-sky/inputs/garand-atmos-1.nc rrtmgp-allsky.nc
	$(RUN_CMD) ./rrtmgp_allsky rrtmgp-allsky.nc ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc 128
	$(RUN_CMD) ./rrtmgp_allsky rrtmgp-allsky.nc ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc 128

check:
	python ./compare-to-reference.py

clean:
	-rm rrtmgp_allsky *.o *.optrpt ../*.optrpt *.mod
