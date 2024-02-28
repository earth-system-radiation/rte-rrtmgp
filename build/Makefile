#!/usr/bin/env make

RTE_DIR           = ../rte-frontend
GAS_OPTICS_DIR    = ../gas-optics
RRTMGP_DIR        = ../rrtmgp-frontend
RTE_KERNEL_DIR    = ../rte-kernels
RRTMGP_KERNEL_DIR = ../rrtmgp-kernels
#
# Compiler variables FC, FCFLAGS must be set in the environment
#
# 
# Make all the libraries though we'll only use the interface + kernels
all:           librte.a                  librrtmgp.a \
               librtekernels.a librtef.a librrtmgpkernels.a librrtmgpf.a
separate-libs: librtekernels.a librtef.a librrtmgpkernels.a librrtmgpf.a

COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<

include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends
include $(RTE_KERNEL_DIR)/Make.depends
include $(RRTMGP_KERNEL_DIR)/Make.depends

#
# If using OpenACC/OpenMP files in *-kernels/accel take precendence  
#
ifeq ($(RTE_KERNELS), accel)
	VPATH = $(RTE_KERNEL_DIR)/accel:$(RRTMGP_KERNEL_DIR)/accel
endif
#
# If using external libraries just compile the interfaces 
#    (and, for now, the OpenACC kernels)
#
ifeq ($(RTE_KERNELS), extern)
	VPATH = $(RTE_KERNEL_DIR)/extern:$(RRTMGP_KERNEL_DIR)/extern:$(RTE_KERNEL_DIR)/accel:$(RRTMGP_KERNEL_DIR)/accel
endif

VPATH +=$(RTE_DIR):$(GAS_OPTICS_DIR):$(RRTMGP_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_KERNEL_DIR)

#
# Library containing Fortran front-end and kernels 
#
librte.a: $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)
	ar -rvs librte.a $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)

#
# Library with just the kernels...
#
librtekernels.a: $(RTE_FORTRAN_KERNELS)
	ar -rvs librtekernels.a $(RTE_FORTRAN_KERNELS)

#
# ... and just the Fortran front-end 
#
librtef.a: $(RTE_FORTRAN_INTERFACE)
	ar -rvs librtef.a $(RTE_FORTRAN_INTERFACE)

#
# As with RTE, libraries with Fortran front-end and kernels, separate and combined 
#
librrtmgp.a: $(RRTMGP_FORTRAN_KERNELS)  $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librrtmgp.a $(RRTMGP_FORTRAN_KERNELS) $(RRTMGP_FORTRAN_INTERFACE)

librrtmgpkernels.a: $(RRTMGP_FORTRAN_KERNELS)
	ar -rvs librrtmgpkernels.a $(RRTMGP_FORTRAN_KERNELS)

librrtmgpf.a: $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librrtmgpf.a $(RRTMGP_FORTRAN_INTERFACE)

clean:
	rm -f *.optrpt *.mod *.o lib*.a
