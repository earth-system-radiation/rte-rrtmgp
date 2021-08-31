#!/bin/bash

unset ARCH
unset CUDA_ARCH
unset CUBHOME

export NCFLAGS="`nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export YAKLHOME="/home/$USER/YAKL"
