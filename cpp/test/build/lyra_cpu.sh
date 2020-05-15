#!/bin/bash

source $MODULESHOME/init/bash
module load gcc hip cmake

unset ARCH
unset CUDA_ARCH
unset CUBHOME
unset ROCPRIMHOME
unset HIPCUBHOME

export NCHOME=/ccs/home/imn/libs_gnu8.1.0/netcdf-4.3.3.1_gnu
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export YAKLHOME="/ccs/home/$USER/YAKL"
