#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load cmake netcdf

export NCHOME=${NETCDF_PATH}
export CC=mpicc
export CXX=mpic++
export CXXFLAGS="-O3"
export ARCH="CUDA"
export CUDA_ARCH="-arch sm_50 --std=c++14 --use_fast_math -O3"
export CUBHOME="/home/$USER/cub"
export YAKLHOME="/home/$USER/YAKL"
