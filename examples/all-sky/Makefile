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

# netcdf Fortran module files has to be in the search path or added via environment variable FCINCLUDE e.g. 
#FCINCLUDE += -I$(NFHOME)/include

# netcdf C and Fortran libraries have to be in the search path or added via environment variable LDFLAGS e.g. 
#LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
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

tests: rrtmgp_allsky
	$(RUN_CMD) bash all_tests.sh  

check:
	$${PYTHON-python} ${RRTMGP_ROOT}/examples/compare-to-reference.py --ref_dir ${RRTMGP_DATA}/examples/all-sky/reference --tst_dir ${RRTMGP_ROOT}/examples/all-sky \
	       --var lw_flux_up lw_flux_dn sw_flux_up sw_flux_dn sw_flux_dir  \
	       --file rrtmgp-allsky-lw.nc rrtmgp-allsky-sw.nc
	$${PYTHON-python} ${RRTMGP_ROOT}/examples/compare-to-reference.py --ref_dir ${RRTMGP_DATA}/examples/all-sky/reference --tst_dir ${RRTMGP_ROOT}/examples/all-sky \
	       --var lw_flux_up lw_flux_dn sw_flux_up sw_flux_dn sw_flux_dir  \
	       --file rrtmgp-allsky-lw-no-aerosols.nc rrtmgp-allsky-sw-no-aerosols.nc

clean:
	-rm rrtmgp_allsky *.o *.optrpt ../*.optrpt *.mod *.nc
