#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc netcdf-c netcdf-fortran cmake cuda/11.5.2 nco git

unset CXXFLAGS

export YAKL_ARCH="CUDA"
export CC=gcc
export CXX=g++
export FC=gfortran
export YAKL_CUDA_FLAGS="-O3 --use_fast_math -DYAKL_AUTO_PROFILE -arch sm_70 -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-O3 -I`nf-config --includedir`"
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKLHOME="/ccs/home/$USER/YAKL"

