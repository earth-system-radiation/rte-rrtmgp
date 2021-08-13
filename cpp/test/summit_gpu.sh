#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc cmake cuda netcdf-cxx4 netcdf nco

export NCINCLUDE="`ncxx4-config --includedir`;`nc-config --includedir`" # must be semi-colon seperated
export NCFLAGS="`ncxx4-config --libs` `nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export ARCH="CUDA"
export YAKL_CUDA_FLAGS="-arch sm_70 --use_fast_math -O3"
export CUBHOME="/ccs/home/$USER/cub"
export YAKLHOME="/ccs/home/$USER/YAKL"

# Set to jsrun, nvprof, or empty depending on how you want to run the tests
export RUNCMD="jsrun -n 1 -a 1 -c 1 -g 1"
