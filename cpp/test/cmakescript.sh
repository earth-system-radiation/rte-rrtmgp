#!/bin/bash

# Clean previous build
./cmakeclean.sh

# Configure new build
printf "NetCDF Flags: $NCFLAGS\n\n"
printf "CXXFLAGS: $CXXFLAGS\n\n"
mkdir -p build && cd build
cmake                                  \
  -DCMAKE_CXX_FLAGS="$CXXFLAGS"        \
  -DNCFLAGS="$NCFLAGS"                 \
  -DARCH="$ARCH"                       \
  -DCUDA_FLAGS="$CUDA_ARCH"            \
  -DYAKL_CUB_HOME="$CUBHOME"           \
  -DYAKL_HIPCUB_HOME="$HIPCUBHOME"     \
  -DYAKL_ROCPRIM_HOME="$ROCPRIMHOME"   \
  -DYAKL_HOME="$YAKLHOME"              \
  ..

