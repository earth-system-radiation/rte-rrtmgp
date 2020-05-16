#!/bin/bash

source $MODULESHOME/init/bash
module load gcc hip cmake

unset ARCH
unset CUDA_ARCH
unset CUBHOME
unset ROCPRIMHOME
unset HIPCUBHOME

export NCFLAGS="`/ccs/home/${USER}/libs_gnu8.1.0/netcdf-4.3.3.1_gnu/bin/ncxx4-config --libs` `/ccs/home/${USER}/libs_gnu8.1.0/netcdf-4.3.3.1_gnu/bin/nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3 -I/ccs/home/${USER}/libs_gnu8.1.0/netcdf-4.3.3.1_gnu/include"
export YAKLHOME="/ccs/home/$USER/YAKL"
