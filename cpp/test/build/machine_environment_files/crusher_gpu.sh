#!/bin/bash

source $MODULESHOME/init/bash
module load PrgEnv-amd cray-hdf5 cray-netcdf cmake craype-accel-amd-gfx90a

unset CXXFLAGS

export YAKL_ARCH="HIP"
export CC=gcc
export CXX=hipcc
export FC=gfortran
export YAKL_HIP_FLAGS="-DYAKL_PROFILE -O3 -ffast-math -DHIP_FAST_MATH -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a -x hip  -I`nc-config --includedir`"
export YAKL_F90_FLAGS="-O3 -I/opt/cray/pe/netcdf/4.8.1.3/gnu/9.1/include"
export CXX_LINK="-L/opt/cray/pe/netcdf/4.8.1.3/amd/4.3/lib -lnetcdf --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib -lamdhip64"
export F90_LINK="-L/opt/cray/pe/netcdf/4.8.1.3/amd/4.3/lib -lnetcdff -lnetcdf"
export YAKLHOME="/ccs/home/$USER/YAKL"

