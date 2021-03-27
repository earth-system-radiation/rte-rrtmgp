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
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

#
# General Purpose Timing Library https://jmrosinski.github.io/GPTL/
# Set environment variable GPTL_DIR to the root of a GPTL installation to build
#   the RFMIP example with timers
#
ifneq ($(origin GPTL_DIR),undefined)
	#
	# Timing library
	#
	FCINCLUDE += -I$(GPTL_DIR)/include
	# Compiler specific
	FCFLAGS   += -DUSE_TIMING
	LDFLAGS   += -L$(GPTL_DIR)/lib
	LIBS      += -lgptl
endif

VPATH = ../

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

#
# Ancillary codes
#
ADDITIONS = mo_simple_netcdf.o mo_rfmip_io.o mo_load_coefficients.o

all: rrtmgp_rfmip_lw rrtmgp_rfmip_sw

rrtmgp_rfmip_lw:   rrtmgp_rfmip_lw.o   $(ADDITIONS) $(RRTMGP_BUILD)/librte.a $(RRTMGP_BUILD)/librrtmgp.a

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

rrtmgp_rfmip_sw:   rrtmgp_rfmip_sw.o   $(ADDITIONS) $(RRTMGP_BUILD)/librrtmgp.a $(RRTMGP_BUILD)/librte.a

rrtmgp_rfmip_sw.o: rrtmgp_rfmip_sw.F90 $(ADDITIONS)

mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o

tests: multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc # Also the template files
	# Files need to have been generated/downloaded before
	$(RUN_CMD) ./rrtmgp_rfmip_lw 8 multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_ROOT}/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
	$(RUN_CMD) ./rrtmgp_rfmip_sw 8 multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc ${RRTMGP_ROOT}/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc

check:
	cd ${RRTMGP_ROOT}/examples/rfmip-clear-sky && python ./compare-to-reference.py --fail=7.e-4

multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc:
	python stage_files.py

clean:
	-rm rrtmgp_rfmip_sw rrtmgp_rfmip_lw *.o *.mod *.optrpt
