#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc/11.2.0 netcdf-c netcdf-fortran cmake nco git

unset CXXFLAGS
unset YAKL_ARCH

export CC=gcc
export CXX=g++
export FC=gfortran
export YAKL_CXX_FLAGS="-O3 -DYAKL_PROFILE -DRRTMGP_CPU_KERNELS -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-O3 -I`nf-config --includedir`"
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKLHOME="/ccs/home/$USER/YAKL"

