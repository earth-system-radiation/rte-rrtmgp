---
layout: "page"
title: "How-to guides"
---
# How-to guides will live here

## How-to: build, run, and test the libraries, examples, and unit-testing codes.

1. Set environment variables `FC` (the Fortran 2003 compiler) and `FCFLAGS` (compiler flags). Examples are provided in the `Compiler-flags.md` file.
2. Set environment variables `RRTMGP_ROOT` to the top-level RTE+RRTMGP directory and `RTE_KERNELS` to `accel` if you want the OpenACC/OpenMP kernels rather than the default.
3. `make libs` in the top-level directory will make the RTE and RRTMGP libraries.
4. The examples and testing codes use netCDF. Set the variables `NCHOME` and `NFHOME` to the roots of the C and Fortran netCDF installations. 
5. Download the RRTMGP data either by cloning the [data repository](https://github.com/earth-system-radiation/rrtmgp-data) or from the [Zenodo archive](https://doi.org/10.5281/zenodo.7988260). Set the environment variable `RRTMGP_DATA` to the root of this directory. 
6. `make tests` to will build and run the test. 
7. Evaluating the results of the tests requires `Python` and the packages described in `environment.yml`. Comparisons can be made with `make check` in the top level directory. 
8. `make` invoked without a target in the top level attempts all three steps.



### Building and testing using (Gnu) make

Sergey Kosukhin and his colleagues at the Max Planck Institute for Meteorology
maintain the `autoconf` branch which adds Gnu `autotools` building to `main` branch.
