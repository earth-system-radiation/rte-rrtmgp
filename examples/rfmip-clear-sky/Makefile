# Location of RTE+RRTMGP libraries, module files.
RRTMGP_DIR = ../../build
# Sets macros FC, FCFLAGS consistent with RTE+RRTMGP
include $(RRTMGP_DIR)/Makefile.conf
#
# RRTMGP library, module files
#
LDFLAGS   += -L$(RRTMGP_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(RRTMGP_DIR)

#
# netcdf library, module files
# C and Fortran interfaces respectively 
#
NCHOME = /opt/local
NFHOME = /Users/robert/Applications/$(FC)
FCINCLUDE += -I$(NFHOME)/include
LDFLAGS   += -L$(NFHOME)/lib -L$(NCHOME)/lib
LIBS      += -lnetcdff -lnetcdf

# Compilation rules
%.o: %.F90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)


#
# Ancillary codes
#
ADDITIONS = mo_rfmip_io.o mo_load_coefficients.o

all: rrtmgp_rfmip_lw

rrtmgp_rfmip_lw: rrtmgp_rfmip_lw.o $(ADDITIONS)

rrtmgp_rfmip_lw.o: rrtmgp_rfmip_lw.F90 $(ADDITIONS)

clean:
	-rm rrtmgp_rfmip_lw *.o *.mod *.optrpt ../*.optrpt
