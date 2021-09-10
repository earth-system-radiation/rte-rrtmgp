#!/bin/bash

export ARCH="CUDA"
unset CXXFLAGS

export NCINCLUDE="`/opt/netcdf_gnu/bin/nc-config --includedir`"
export NCFLAGS="`/opt/netcdf_gnu/bin/nc-config --libs`"
export CC=gcc
export CXX=g++
export YAKL_CUDA_FLAGS="-arch sm_35 --std=c++14 --use_fast_math -O3"
export YAKLHOME="/home/$USER/YAKL"
