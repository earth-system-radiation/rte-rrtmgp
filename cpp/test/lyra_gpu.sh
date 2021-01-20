#!/bin/bash

source $MODULESHOME/init/bash
module load gcc hip cmake

export NCFLAGS="`/ccs/home/imn/libs_hip/netcdf-4.3.3.1/bin/ncxx4-config --libs` `/ccs/home/imn/libs_gnu8.1.0/netcdf-4.3.3.1_gnu/bin/nc-config --libs`"
export CC=gcc
export CXX=hipcc
export CXXFLAGS="-O3 -I/ccs/home/imn/libs_hip/netcdf-4.3.3.1/include -I/ccs/home/imn/libs_gnu8.1.0/netcdf-4.3.3.1_gnu/include"
export ARCH="HIP"
export YAKLHOME="/ccs/home/$USER/YAKL"
export HIPCUBHOME="/ccs/home/$USER/hipCUB"
export ROCPRIMHOME="/ccs/home/$USER/rocPRIM"
