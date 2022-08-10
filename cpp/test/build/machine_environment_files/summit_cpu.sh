#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc cmake cuda netcdf-c nco

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCINCLUDE="`nc-config --includedir`"
export NCFLAGS="`nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export YAKL_CUDA_FLAGS=""
export YAKLHOME="/ccs/home/$USER/YAKL"

# Set to jsrun, nvprof, or empty depending on how you want to run the tests
export RUNCMD="jsrun -n 1 -a 1 -c 1 -g 1"
