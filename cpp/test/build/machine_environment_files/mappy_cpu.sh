unset ARCH
unset CUDA_ARCH
unset CUBHOME

# Just use the scream env

this_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

export CC=gcc
export CXX=g++
export FC=gfortran
export CXX_LINK="`nc-config --libs`"
export F90_LINK="`nf-config --flibs`"
export YAKLHOME=$this_dir/../../../../../../../../../../externals/YAKL
export KOKKOSHOME=$this_dir/../../../../../../../../../../externals/ekat/extern/kokkos
