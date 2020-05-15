#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load cmake netcdf

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCHOME=${NETCDF_PATH}
export CC=mpicc
export CXX=mpic++
export CXXFLAGS="-O3"
export YAKLHOME="/home/$USER/YAKL"
