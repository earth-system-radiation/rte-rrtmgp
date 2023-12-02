---
layout: "page"
title: "How to build the libraries, tests, and examples, run the tests, and verify the results"
---
## In a nutshell 
In the root directory: 
- `make libs`  makes the RTE and RRTMGP libraries, the unit tests, and the examples 
- `make tests` runs the tests
- `make check` uses Python to verify results against reference calculations
- `make` invoked without a target in the top level attempts all three steps.

Evaluating the results of the tests requires `Python` and the packages described in `environment*.yml`.

## Building and testing using the handbuilt Makefiles 

Before using the Makefiles supplied with the `RTE+RRTMGP` repository, the environment variables `FC` and
`FCFLAGS`, identifying the Fortran compiler and flags passed to it, need to be set. 

To build any of the examples in `examples/` or `tests` the locations of the C and Fortran netCDF libraries and the 
location of the netCDF Fortran module file (`netcdf.mod`) must be in the search path. 
Non-standard paths can also be added via macros `FCINCLUDE` and/or `LDFLAGS`.

## Building and testing using (Gnu) autotools 

Sergey Kosukhin and his colleagues at the Max Planck Institute for Meteorology
maintain the `autoconf` branch which adds Gnu `autotools` building to `main` branch.

## Supplying data 

Running the tests and verifying the results requires the RRTMGP data. Clone the 
[data repository](https://github.com/earth-system-radiation/rrtmgp-data) or download the 
[Zenodo archive](https://doi.org/10.5281/zenodo.7988260). Set the environment variable `RRTMGP_DATA` 
to the root of this directory. 

## Example compiler flags 

### Gnu Fortran 
(see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/continuous-integration.yml))
`FC: `gfortran-10` or `gfortran-11` or `gfortran-12`
#### Debugging flags
`FCFLAGS: "-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -finit-real=nan -DRTE_USE_CBOOL"`  
#### Even stricter debugging flags
`FCFLAGS: "-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -fbacktrace -finit-real=nan -DRTE_USE_CBOOL -pedantic -g -Wall"`  

### Intel Fortran Classic 
(see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/containerized-ci.yml))
`FC: ifort`  
#### Debugging flags
`FCFLAGS: "-m64 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132 -check bounds,uninit,pointers,stack -stand f08"`  
#### Optimization flags:  
`FCFLAGS:"-m64 -O3 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132"`

### Intel Fortran 
(LLVM, see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/containerized-ci.yml))
`FC: ifort`  
#### Debugging flags
`FCFLAGS: "-debug -traceback -heap-arrays -assume realloc_lhs -extend-source 132 -stand f08"`  
#### Using OpenMP GPU offload 
See [this open issue](https://github.com/earth-system-radiation/rte-rrtmgp/issues/194)

### NVFortran
(see also the see also the [continuous integration](https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/.github/workflows/containerized-ci.yml))
`FC: nvfortran`
#### Debugging flags
`FCFLAGS: "-g -Minfo -Mbounds -Mchkptr -Mstandard -Kieee -Mchkstk -Mallocatable=03  -Mpreprocess"`
#### Optimization flags:  
`FCFLAGS: "-O3 -fast -Minfo -Mallocatable=03 -Mpreprocess"`

### HPE CCE for GPU using OpenMP-acc: crayftn   -- requires at least CCE 14.0.0
`FC: crayftn`
#### Debugging flags  (these appear to be insufficient during the link stage)
`FCFLAGS: "-hnoacc -homp -O0"`