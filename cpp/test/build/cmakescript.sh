#!/bin/bash

savedir=`pwd`
cd ${YAKLHOME}
git rev-parse HEAD >& $savedir/../../yakl-git-hash.txt
cd $savedir

this_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

# Clean previous build
$this_dir/cmakeclean.sh

# Configure new build
cmake                                          \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS}"   \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS}"       \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS}"         \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS}"       \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"         \
  -DYAKL_ARCH="$YAKL_ARCH"                     \
  -DYAKL_HOME="$YAKLHOME"                      \
  -DKokkos_DIR="$KOKKOSHOME"                   \
  -DCXX_LINK="$CXX_LINK"                       \
  -DF90_LINK="$F90_LINK"                       \
  -DCMAKE_Fortran_MODULE_DIRECTORY="`pwd`/fortran_module_files" \
  -DCMAKE_CUDA_ARCHITECTURES=OFF               \
  $this_dir/..
