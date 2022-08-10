#!/bin/bash

module load/9.2.0
module load cuda11.1/toolkit/11.1.1
spack load netcdf-cxx4%gcc

export NCINCLUDE="`nc-config --includedir`"
export NCFLAGS="`nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
unset ARCH
unset CUDA_ARCH
unset CUBHOME
export YAKLHOME="/home/users/$USER/YAKL"
