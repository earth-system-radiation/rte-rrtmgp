#!/usr/bin/env make

RTE_DIR = ../rte
RRTMGP_DIR = ../rrtmgp
#
# Choose kernels depending on environment variable
#
ifeq ($(RTE_KERNELS), openacc)
  RTE_KERNEL_DIR    = ../rte/kernels-openacc
  RRTMGP_KERNEL_DIR = ../rrtmgp/kernels-openacc
endif
RTE_KERNEL_DIR += ../rte/kernels
RRTMGP_KERNEL_DIR += ../rrtmgp/kernels
VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)
#
# Compiler variables FC, FCFLAGS must be set in the environment
#
# Make all the libraries though we'll only use the interface + kernels
all:           librte.a                        librrtmgp.a
separate-libs: librtekernels.a librteinfacef.a librrtmgp.a

COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<

include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends

librte.a: $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)
	ar -rvs librte.a $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)

librtekernels.a: $(RTE_FORTRAN_KERNELS)
	ar -rvs librtekernels.a $(RTE_FORTRAN_KERNELS)

librteinfacef.a: $(RTE_FORTRAN_INTERFACE)
	ar -rvs librteinfacef.a $(RTE_FORTRAN_INTERFACE)

librrtmgp.a: $(RRTMGP_SRC)
	ar -rvs librrtmgp.a $(RRTMGP_SRC)

clean:
	rm -f *.optrpt *.mod *.o librrtmgp.a librte.a
