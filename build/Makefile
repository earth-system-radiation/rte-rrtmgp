#!/bin/make
RTE_DIR = ../rte
RTE_KERNEL_DIR = ../rte/kernels
RRTMGP_DIR = ../rrtmgp
RRTMGP_KERNEL_DIR = ../rrtmgp/kernels
VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)

#
# Setting macro TIMING=yes uses routines from the General Purpose Timing Library
#  https://jmrosinski.github.io/GPTL/
#
TIMING=yes
# Compiler specific - e.g. turn off for pgfortran, set to -cpp for gfortran
FCFLAGS += -fpp
ifeq ($(TIMING),yes)
	#
	# Timing library
	#
	FCINCLUDE += -I$(TIME_DIR)/include
	# Compiler specific
	FCFLAGS += -DUSE_TIMING
	# Installation specific
	TIME_DIR = $(HOME)/Codes/GPTL/gptl-v5.5.3/macos
	LDFLAGS   += -L$(TIME_DIR)/lib
	LIBS      += -lgptl
endif

include Makefile.conf
include Makefile.rules

all: librrtmgp.a librte.a

include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends

librte.a: $(RTE_SRC)
	ar -rvs librte.a $(RTE_SRC)

librrtmgp.a: $(RRTMGP_SRC)
	ar -rvs librrtmgp.a $(RRTMGP_SRC)

clean:
	rm -f *.optrpt *.mod *.o librrtmgp.a librte.a
