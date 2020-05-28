# Compiler flag Examples

Before using the Makefiles supplied with the `RTE+RRTMGP` repository, the environment variables `FC` and
`FCFLAGS`, identifying the Fortran compiler and flags passed to it, need to be set. Here are some examples
used during development and testing.

To build and of the executables in `examples/` or `tests` the locations of the C and Fortran netCDF libraries
need to be set via environment variables `NCHOME` and `NFHOME`.

## Gnu Fortran
`FC: gfortran-8` or `gfortran-9`  
Debugging flags:  
`FFLAGS: "-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -finit-real=nan -DUSE_CBOOL"`  
Optimization flags:  
`FCFLAGS: "-ffree-line-length-none -m64 -std=f2008 -march=native -pedantic -g -fbounds-check -Wall -fbacktrace -finit-real=nan -DUSE_CBOOL"`  

## Intel Fortran
`FC: ifort`  
Debugging flags:  
`FFLAGS: "-m64 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132 -check bounds,uninit,pointers,stack -stand f08"`  
Optimization flags  
`FCFLAGS:"-m64 -O3 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132"`

## PGI Fortran
`FC: pgfortran`

Debugging flags:

`FFLAGS: "-g -Minfo -Mbounds -Mchkptr -Mstandard -Kieee -Mchkstk -Mallocatable=03  -Mpreprocess"`

Optimization flags:

`FCFLAGS: "-m64 -O3 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132"`
