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

#
# Setting macro TIMING=yes uses routines from the General Purpose Timing Library
#  https://jmrosinski.github.io/GPTL/
# If TIMING=yes evironment variables TIME_DIR must be set to root of GPTL installation
#
TIMING=no
# Compiler specific - e.g. turn off for pgfortran, set to -cpp for gfortran
#FCFLAGS += -fpp
ifeq ($(TIMING),yes)
	#
	# Timing library
	#
	FCINCLUDE += -I$(TIME_DIR)/include
	# Compiler specific
	FCFLAGS += -DUSE_TIMING
	LDFLAGS   += -L$(TIME_DIR)/lib
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

rrtmgp_rfmip_sw:   rrtmgp_rfmip_sw.o   $(ADDITIONS) $(RRTMGP_BUILD)/librte.a $(RRTMGP_BUILD)/librrtmgp.a

rrtmgp_rfmip_sw.o: rrtmgp_rfmip_sw.F90 $(ADDITIONS)

mo_rfmip_io.o:          mo_rfmip_io.F90          mo_simple_netcdf.o

mo_load_coefficients.o: mo_load_coefficients.F90 mo_simple_netcdf.o

clean:
	-rm rrtmgp_rfmip_sw rrtmgp_rfmip_lw *.o *.mod *.optrpt
