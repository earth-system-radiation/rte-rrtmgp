#!/bin/bash

unset CXXFLAGS

export YAKL_ARCH="CUDA"
export CC=gcc
export CXX=g++
export YAKL_CUDA_FLAGS="-O3 --use_fast_math -arch sm_50 -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-O3 -I`nf-config --includedir`"
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKLHOME="/home/$USER/YAKL"

