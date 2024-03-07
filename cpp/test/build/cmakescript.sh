#!/bin/bash -ex

savedir=`pwd`
cd ${YAKLHOME}
git rev-parse HEAD >& $savedir/../../yakl-git-hash.txt
cd $savedir

this_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

# Clean previous build
$this_dir/cmakeclean.sh

# Configure new build. This will pass along arguments to this script
# to CMake so you can set things like CMAKE_BUILD_TYPE easily.
#
# Example: Enable Kokkos and debug build type:
#  $path_to_build/cmakescript.sh -DRRTMGP_ENABLE_KOKKOS=On -DCMAKE_BUILD_TYPE=Debug
#
# Example: Enable Kokkos and YAKL with debug build type:
#  $path_to_build/cmakescript.sh -DRRTMGP_ENABLE_KOKKOS=On -DRRTMGP_ENABLE_YAKL=On -DCMAKE_BUILD_TYPE=Debug
#
# Example: Enable Kokkos with release build type:
#  $path_to_build/cmakescript.sh -DRRTMGP_ENABLE_KOKKOS=On -DCMAKE_BUILD_TYPE=Release
#
# Example: Enable Kokkos with release build type and machine setup
#  $path_to_build/cmakescript.sh -C $KOKKOSHOME/../../cmake/machine-files/mappy.cmake -DRRTMGP_ENABLE_KOKKOS=On -DCMAKE_BUILD_TYPE=Release
#
# Example: Enable YAKL with release build type with threads
#  YAKL_ARCH=OPENMP $path_to_build/cmakescript.sh -DRRTMGP_ENABLE_YAKL=On -DCMAKE_BUILD_TYPE=Release
cmake                                          \
  $@                                           \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS}"   \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS}"       \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS}"         \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS}"       \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"         \
  -DYAKL_ARCH="$YAKL_ARCH"                     \
  -DYAKL_HOME="$YAKLHOME"                      \
  -DKokkos_DIR="$KOKKOSHOME"                   \
  ${KOKKOS_CONFIG}                             \
  -DCXX_LINK="$CXX_LINK"                       \
  -DF90_LINK="$F90_LINK"                       \
  -DCMAKE_Fortran_MODULE_DIRECTORY="`pwd`/fortran_module_files" \
  -DCMAKE_CUDA_ARCHITECTURES=OFF               \
  $this_dir/..
