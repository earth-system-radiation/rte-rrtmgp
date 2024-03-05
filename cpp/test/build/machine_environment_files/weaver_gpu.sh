#!/bin/bash

unset ARCH
unset CUDA_ARCH
unset CUBHOME

# Just use the scream env

this_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

export YAKL_ARCH="CUDA"
export CC=gcc
export CXX=g++
export FC=gfortran
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKL_CUDA_FLAGS="-arch sm_70 -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-I`nf-config --includedir`"
E3SM_ROOT=$this_dir/../../../../../../../../../..
export YAKLHOME=$E3SM_ROOT/externals/YAKL
export KOKKOSHOME=$E3SM_ROOT/externals/ekat/extern/kokkos
export KOKKOS_CONFIG="-C $E3SM_ROOT/externals/ekat/cmake/machine-files/weaver.cmake"
