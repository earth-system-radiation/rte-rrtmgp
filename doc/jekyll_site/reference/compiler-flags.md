---
layout: page
title: Example compiler flags
---

Here are sets of compiler flags used to compile and test the code over time. They are provided as a starting point - other configurations may work as well or better.

In these examples `FC` is the Fortran compilers using flags `FCFLAGS`

# Gnu Fortran

(see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/continuous-integration.yml))
`FC`: `gfortran-12` or `gfortran-13` or `gfortran-14`

## Debugging flags

`FCFLAGS: "-ffree-line-length-none -std=f2008 -fbounds-check -finit-real=nan"`

## Even stricter debugging flags

`FCFLAGS: "-ffree-line-length-none -std=f2008 -fbounds-check -finit-real=nan -fbacktrace -pedantic -g -Wall"`

# Intel Fortran Classic

(see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/containerized-ci.yml))
`FC: ifort`

## Debugging flags

`FCFLAGS: "-m64 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132 -check bounds,uninit,pointers,stack -stand f08"`

## Optimization flags:

`FCFLAGS:"-m64 -O3 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132"`

# Intel Fortran

(LLVM, see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/containerized-ci.yml))
`FC: ifx`

## Debugging flags

`FCFLAGS: "-debug -traceback -heap-arrays -assume realloc_lhs -extend-source 132 -stand f08"`

## Using OpenMP GPU offload

See [this open issue](https://github.com/earth-system-radiation/rte-rrtmgp/issues/194)

# NVFortran

(see also the see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/containerized-ci.yml))
`FC: nvfortran`

## Debugging flags

`FCFLAGS: "-g -Minfo -Mbounds -Mchkptr -Mstandard -Kieee -Mchkstk -Mallocatable=03  -Mpreprocess"`

## Optimization flags:

`FCFLAGS: "-O3 -fast -Minfo -Mallocatable=03 -Mpreprocess"`

# HPE CCE for GPU using OpenMP-acc: crayftn -- requires at least CCE 14.0.0

`FC: crayftn`

## Debugging flags (these appear to be insufficient during the link stage)

`FCFLAGS: "-hnoacc -homp -O0"`
