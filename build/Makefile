#!/usr/bin/env make

RTE_DIR           = ../rte-frontend
GAS_OPTICS_DIR    = ../gas-optics
RRTMGP_DIR        = ../rrtmgp-frontend
RTE_KERNEL_DIR    = ../rte-kernels
RRTMGP_KERNEL_DIR = ../rrtmgp-kernels
#
# Compiler variables FC, FCFLAGS must be set in the environment
#
# Make all the libraries though we'll only use the interface + kernels
all:           librte.a                  librrtmgp.a		librtecudakernels.a librrtmgpcudakernels.a
separate-libs: librtekernels.a librtef.a librrtmgpkernels.a librrtmgpf.a

COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
%.o: %.F90
	$(COMPILE) $<

COMPILE_CU = nvc++ -c -std=c++17 -O3 -DRTE_RRTMGP_GPU_MEMPOOL_CUDA
%.o: %.cu
	$(COMPILE_CU) $<

include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends
include $(RTE_KERNEL_DIR)/Make.depends
include $(RRTMGP_KERNEL_DIR)/Make.depends

include ../cuda/Make.depends

VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR):$(GAS_OPTICS_DIR)
#
# If using OpenACC/OpenMP files in *-kernels/accel take precendence  
#
ifeq ($(RTE_KERNELS), accel)
	VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR)/accel:$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)/accel:$(RRTMGP_KERNEL_DIR):$(GAS_OPTICS_DIR)
endif
ifeq ($(RTE_KERNELS), extern)
	VPATH = ../cuda:$(RTE_DIR):$(RTE_KERNEL_DIR)/accel:$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)/extern:$(RRTMGP_KERNEL_DIR)/accel:$(RRTMGP_KERNEL_DIR):$(GAS_OPTICS_DIR)
endif

librte.a: $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)
	ar -rvs librte.a $(RTE_FORTRAN_KERNELS) $(RTE_FORTRAN_INTERFACE)

librtecudakernels.a: $(RTE_CUDA_KERNELS)
	ar -rvs librtecudakernels.a $(RTE_CUDA_KERNELS)

librtef.a: $(RTE_FORTRAN_INTERFACE)
	ar -rvs librtef.a $(RTE_FORTRAN_INTERFACE)

librrtmgp.a: $(RRTMGP_CUDA_KERNELS) $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librrtmgp.a $(RRTMGP_CUDA_KERNELS) $(RRTMGP_FORTRAN_INTERFACE)

librrtmgpcudakernels.a: $(RRTMGP_CUDA_KERNELS)
	ar -rvs librrtmgpcudakernels.a $(RRTMGP_CUDA_KERNELS)

librrtmgpkernels.a: $(RRTMGP_FORTRAN_KERNELS)
	ar -rvs librrtmgpkernels.a $(RRTMGP_FORTRAN_KERNELS)

librrtmgpf.a: $(RRTMGP_FORTRAN_INTERFACE)
	ar -rvs librrtmgpf.a $(RRTMGP_FORTRAN_INTERFACE)


clean:
	rm -f *.optrpt *.mod *.o lib*.a
