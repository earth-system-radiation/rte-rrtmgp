#!/bin/bash

unset CXXFLAGS

export YAKL_ARCH="CUDA"
export CC=gcc-11
export CXX=g++-11
export YAKL_CUDA_FLAGS="-O3 -arch sm_75 -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-O3 -I`nf-config --includedir`"
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKLHOME="/home/$USER/YAKL"

