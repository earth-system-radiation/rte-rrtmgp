#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load DefApps gcc cmake cuda netcdf-cxx4 netcdf

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCFLAGS="`ncxx4-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export YAKLHOME="/ccs/home/$USER/YAKL"
