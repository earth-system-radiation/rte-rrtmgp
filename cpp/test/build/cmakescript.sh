#!/bin/bash

# Clean previous build
./cmakeclean.sh

# Configure new build
cmake                                    \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}"   \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS}" \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"   \
  -DYAKL_ARCH="$YAKL_ARCH"               \
  -DYAKL_HOME="$YAKLHOME"                \
  -DCXX_LINK="$CXX_LINK"                 \
  -DF90_LINK="$F90_LINK"                 \
  -DCMAKE_Fortran_MODULE_DIRECTORY="`pwd`/fortran_module_files" \
  ..
