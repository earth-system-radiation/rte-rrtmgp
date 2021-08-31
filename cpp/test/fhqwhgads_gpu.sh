#!/bin/bash

export NCFLAGS="`nc-config --libs`"
export CC=gcc
export CXX=g++
export CXXFLAGS="-O3"
export ARCH="CUDA"
export CUDA_ARCH="-arch sm_50 --use_fast_math -O3"
export CUBHOME="/home/$USER/cub"
export YAKLHOME="/home/$USER/YAKL"
