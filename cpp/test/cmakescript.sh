#!/bin/bash

# Clean previous build
./cmakeclean.sh

# Configure new build
printf "NetCDF Flags: $NCFLAGS\n\n"
printf "NetCDF Include: $NCINCLUDE\n\n"
printf "CXXFLAGS: $CXXFLAGS\n\n"
mkdir -p build && cd build
cmake                                    \
  -DCMAKE_CXX_FLAGS="$CXXFLAGS"          \
  -DNCINCLUDE="$NCINCLUDE"               \
  -DNCFLAGS="$NCFLAGS"                   \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}"   \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS}" \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}" \
  -DYAKL_ARCH="$ARCH"                    \
  -DYAKL_HOME="$YAKLHOME"                \
  -DFORTRAN_LINK="`nf-config --flibs`"         \
  ..
