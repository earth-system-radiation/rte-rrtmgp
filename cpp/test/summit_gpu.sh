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
export CUDA_ARCH="-arch sm_70 --std=c++14 --use_fast_math -O3"
export YAKLHOME="/ccs/home/$USER/YAKL"
