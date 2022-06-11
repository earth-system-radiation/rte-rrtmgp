---
layout: "page"
title: "How-to guides"
---
# How-to guides will live here

## How to: build and test the libraries and examples

### Building and testing using (Gnu) make

1. Set environment variables `FC` (the Fortran 2003 compiler) and `FCFLAGS` (compiler flags).
   Examples are provided in the `Compiler-flags.md` file.
2. Set environment variable `RRTMGP_ROOT` to the top-level RTE+RRTMGP directory.
   Set the variables `NCHOME` and `NFHOME` to the roots of the C and Fortran
   netCDF installations. (Building the libraries alone )
3. Set environment variable `RTE_KERNELS` to `openacc` if you want the OpenACC/OpenMP
   kernels rather than the default.
4. `make libs` in the top-level directory will make the RTE and RRTMGP libraries
   and the regression tests in in `examples/` and `tests/`. Libraries and module
   files are in `build/`; examples and tests are in the subdirectory containing
   their source code.
5. `make tests` runs the examples and regression tests.
   (A few files need to be downloaded for `examples/rfmip-clear-sky`. The default
   is to download with a Python script is shell script using `wget` is also available.)
6. Comparisons can be made with `make check` in the top level directory.
   Evaluating the results of the tests requires `Python` and the packages
   described in `environment.yml`. One approach is to use
   `conda env create -f environment.yml; conda activate rte_rrtmgp_test; make check`  
7. `make` invoked without a target in the top level attempts all three steps.

### Building and testing using (Gnu) make

Sergey Kosukhin and his colleagues at the Max Planck Institute for Meteorology
maintain the `autoconf` branch which adds Gnu `autotools` building to `main` branch.
