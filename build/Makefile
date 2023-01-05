#!/usr/bin/env make

RTE_DIR = ../rte-frontend
RRTMGP_DIR = ../rrtmgp

RTE_KERNEL_DIR += ../rte-kernels
RRTMGP_KERNEL_DIR += ../rrtmgp/kernels
VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)
#
# Compiler variables FC, FCFLAGS must be set in the environment
#
# Make all the libraries though we'll only use the interface + kernels
all:           librte.a                  librrtmgp.a
separate-libs: librtekernels.a librtef.a librrtmgpkernels.a librrtmgpf.a

COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<

include $(RTE_DIR)/Make.depends
include $(RTE_KERNEL_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends

#
# Choose kernels depending on environment variable
#
ifeq ($(RTE_KERNELS), openacc)
  RTE_KERNEL_DIR    = $(RTE_KERNEL_DIR)/accel:$(RTE_KERNEL_DIR)
  RRTMGP_KERNEL_DIR = ../rrtmgp/kernels-openacc
endif

librte.a: $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)
	ar -rvs librte.a $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)

librtekernels.a: $(RTE_FORTRAN_KERNELS)
	ar -rvs librtekernels.a $(RTE_FORTRAN_KERNELS)

librtef.a: $(RTE_FORTRAN_INTERFACE)
	ar -rvs librtef.a $(RTE_FORTRAN_INTERFACE)

librrtmgp.a: $(RRTMGP_FORTRAN_KERNELS)  $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librrtmgp.a $(RRTMGP_FORTRAN_KERNELS) $(RRTMGP_FORTRAN_INTERFACE)

librrtmgpkernels.a: $(RRTMGP_FORTRAN_KERNELS)
	ar -rvs librrtmgpkernels.a $(RRTMGP_FORTRAN_KERNELS)

librrtmgpf.a: $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librrtmgpf.a $(RRTMGP_FORTRAN_INTERFACE)

clean:
	rm -f *.optrpt *.mod *.o lib*.a
