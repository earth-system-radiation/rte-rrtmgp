#!/bin/bash

module purge
module load sems-env sems-python/3.5.2 sems-gcc/9.2.0 sems-cmake/3.12.2 sems-git/2.10.1 sems-openmpi/4.0.2
module load sems-netcdf
export CTEST_PARALLEL_LEVEL=1

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCFLAGS="`ncxx4-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export YAKLHOME="${HOME}/codes/yakl/YAKL"
