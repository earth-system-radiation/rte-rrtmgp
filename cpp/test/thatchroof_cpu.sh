#!/bin/bash

source $MODULESHOME/init/bash
module purge
module load cmake

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCINCLUDE="`ncxx4-config --includedir`;`nc-config --includedir`" # must be semi-colon seperated
export NCFLAGS="`ncxx4-config --libs` `nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export YAKLHOME="/home/$USER/YAKL"
