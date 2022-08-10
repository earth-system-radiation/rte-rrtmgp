#!/bin/bash

unset YAKL_ARCH
unset CXXFLAGS

export CC=gcc-11
export CXX=g++-11
export YAKL_CXX_FLAGS="-O0 -g -DYAKL_DEBUG -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-O0 -g              -I`nf-config --includedir`"
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKLHOME="/home/$USER/YAKL"

