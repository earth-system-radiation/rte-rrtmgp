#
# Meant to be called as a sub-make; relies on variables set in ../../Makefile
#
COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<
%: %.o
		$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

VPATH = ./:$(RRTMGP_ROOT)/examples/:$(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/extensions/aerosol_optics:$(TEST_ROOT)/test/util/src/

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_aerosol_optics.o mo_test_files_io.o mo_load_aerosol_coefficients.o  mo_simple_netcdf.o

#
# Targets
#
all: test_aerosol_optics

test_aerosol_optics:      $(ADDITIONS) test_aerosol_optics.o

test_aerosol_optics.o: $(ADDITIONS) test_aerosol_optics.F90

mo_load_aerosol_coefficients.o: mo_simple_netcdf.o mo_aerosol_optics.o mo_load_aerosol_coefficients.F90

mo_test_files_io.o: mo_simple_netcdf.o mo_test_files_io.F90

clean:
	-rm ../test_aerosol_optics *.o *.optrpt ../*.optrpt *.mod
