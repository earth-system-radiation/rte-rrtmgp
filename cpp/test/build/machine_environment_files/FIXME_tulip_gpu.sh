#!/bin/bash

module load/9.2.0
module load cuda11.1/toolkit/11.1.1
spack load netcdf-cxx4%gcc

export NCINCLUDE="`nc-config --includedir`"
export NCFLAGS="`nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export ARCH="CUDA"
export CUDA_ARCH="-arch sm_70 --std=c++14 --use_fast_math -O3"
export CUBHOME="/cm/shared/apps/cuda11.0/toolkit/11.0.3/include/cub"
export YAKLHOME="/home/users/$USER/YAKL"
