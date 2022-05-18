#!/bin/bash

unset ARCH
unset CXXFLAGS

export NCINCLUDE="`nc-config --includedir`"
export NCFLAGS="`nc-config --libs`"
export CC=gcc
export CXX=g++
export FC=gfortran
export YAKL_CXX_FLAGS="-O3 -march=native -DYAKL_PROFILE"
export YAKL_F90_FLAGS="-O3 -march=native -I/usr/include"
export YAKLHOME="/home/$USER/YAKL"
