#!/bin/bash

./cmakeclean.sh
printf "NetCDF Flags: $NCFLAGS\n\n"

############################################################################
## RUN THE CONFIGURE
############################################################################
CXXFLAGS="$CXXFLAGS -std=c++14"

printf "CXXFLAGS: $CXXFLAGS\n\n"

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


