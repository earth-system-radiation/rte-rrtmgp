#!/bin/bash

unset ARCH

export NCINCLUDE="`nc-config --includedir`"
export NCFLAGS="`nc-config --libs`"
export CC=gcc-11
export CXX=g++-11
export FC=gfortran-11
export YAKLHOME="/home/$USER/YAKL"

export YAKL_CXX_FLAGS="-O0 -g -DYAKL_DEBUG"
export YAKL_F90_FLAGS="-O0 -g -fcheck=all -I`nf-config --includedir`"

unset CXXFLAGS

