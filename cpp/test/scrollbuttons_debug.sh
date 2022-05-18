#!/bin/bash

unset ARCH
unset CXXFLAGS

export NCINCLUDE="`nc-config --includedir`"
export NCFLAGS="`nc-config --libs`"
export CC=gcc
export CXX=g++
export FC=gfortran
export YAKL_CXX_FLAGS="-O0 -g -DYAKL_DEBUG -DYAKL_PROFILE"
export YAKL_F90_FLAGS="-O0 -g -I/usr/include"
export YAKLHOME="/home/$USER/YAKL"
