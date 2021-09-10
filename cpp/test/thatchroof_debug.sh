#!/bin/bash

unset ARCH
unset CXXFLAGS

export NCINCLUDE="`/opt/netcdf_gnu/bin/nc-config --includedir`"
export NCFLAGS="`/opt/netcdf_gnu/bin/nc-config --libs`"
export CC=gcc
export CXX=g++
export YAKL_CXX_FLAGS="-O0 -g -DYAKL_DEBUG"
export YAKLHOME="/home/$USER/YAKL"
