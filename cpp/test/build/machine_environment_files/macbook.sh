#!/bin/bash

export CTEST_PARALLEL_LEVEL=1

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCFLAGS="`nc-config --libs`"
export NFFLAGS="`nc-config --flibs` -lnetcdf -lnetcdff"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3 -I`nc-config --includedir`"
export F90FLAGS="-O3 -I`nc-config --fflags`"
export YAKLHOME="${HOME}/codes/yakl/YAKL"
