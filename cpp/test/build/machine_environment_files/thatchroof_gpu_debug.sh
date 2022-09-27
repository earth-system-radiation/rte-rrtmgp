#!/bin/bash

unset CXXFLAGS

export YAKL_ARCH="CUDA"
export CC=gcc-11
export CXX=g++-11
export YAKL_CUDA_FLAGS="-O0 -g -DYAKL_DEBUG -DRRTMGP_EXPENSIVE_CHECKS -DRRTMGP_DEBUG -arch sm_75 -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-O0 -g -I`nf-config --includedir`"
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKLHOME="/home/$USER/YAKL"

