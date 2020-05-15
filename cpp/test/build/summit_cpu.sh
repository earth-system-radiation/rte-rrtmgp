#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc cmake cuda netcdf-cxx4 netcdf

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCHOME=${OLCF_NETCDF_CXX4_ROOT}
export CC=mpicc
export CXX=mpic++
export CXXFLAGS="-O3"
export YAKLHOME="/ccs/home/$USER/YAKL"
