#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load cmake

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCFLAGS="`ncxx4-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export YAKLHOME="/home/$USER/YAKL"
