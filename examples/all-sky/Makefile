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

VPATH = ../:$(RRTMGP_ROOT)/rrtmgp-frontend # Needed for cloud_optics and aerosol_optics

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS  = mo_load_coefficients.o mo_simple_netcdf.o 
ADDITIONS += mo_cloud_optics_rrtmgp.o mo_load_cloud_coefficients.o
ADDITIONS += mo_aerosol_optics_rrtmgp_merra.o mo_load_aerosol_coefficients.o

#
# Targets
#
all: rrtmgp_allsky

rrtmgp_allsky: $(ADDITIONS) rrtmgp_allsky.o

rrtmgp_allsky.o: $(ADDITIONS) rrtmgp_allsky.F90

mo_cloud_optics_rrtmgp.o:         mo_cloud_optics_rrtmgp.F90
mo_aerosol_optics_rrtmgp_merra.o: mo_aerosol_optics_rrtmgp_merra.F90
mo_load_coefficients.o:         mo_simple_netcdf.o                                  mo_load_coefficients.F90
mo_load_cloud_coefficients.o:   mo_simple_netcdf.o mo_cloud_optics_rrtmgp.o         mo_load_cloud_coefficients.F90
mo_load_aerosol_coefficients.o: mo_simple_netcdf.o mo_aerosol_optics_rrtmgp_merra.o mo_load_aerosol_coefficients.F90

tests:
	$(RUN_CMD) ./rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-lw.nc 
	$(RUN_CMD) ./rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc ${RRTMGP_DATA}/rrtmgp-aerosols-merra-sw.nc 
	$(RUN_CMD) ./rrtmgp_allsky 24 72 1 rrtmgp-allsky-lw-no-aerosols.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-lw-g256.nc ${RRTMGP_DATA}/rrtmgp-clouds-lw.nc 
	$(RUN_CMD) ./rrtmgp_allsky 24 72 1 rrtmgp-allsky-sw-no-aerosols.nc \
	    ${RRTMGP_DATA}/rrtmgp-gas-sw-g224.nc ${RRTMGP_DATA}/rrtmgp-clouds-sw.nc  

check:
	python ./compare-to-reference.py --allsky_file rrtmgp-allsky-lw.nc
	python ./compare-to-reference.py --allsky_file rrtmgp-allsky-sw.nc
	python ./compare-to-reference.py --allsky_file rrtmgp-allsky-lw-no-aerosols.nc
	python ./compare-to-reference.py --allsky_file rrtmgp-allsky-sw-no-aerosols.nc

clean:
	-rm rrtmgp_allsky *.o *.optrpt ../*.optrpt *.mod *.nc
