#!/bin/bash

unset ARCH

export NCINCLUDE="`nc-config --includedir`"
export NCFLAGS="`nc-config --libs`"
export CC=gcc-11
export CXX=g++-11
export FC=gfortran-11
export YAKLHOME="/home/$USER/YAKL"

export YAKL_CXX_FLAGS="-O3 -march=native -DYAKL_PROFILE"
export YAKL_F90_FLAGS="-O3 -march=native -I`nf-config --includedir`"

unset CXXFLAGS

