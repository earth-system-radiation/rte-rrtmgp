# Compiler flag Examples

Before using the Makefiles supplied with the `RTE+RRTMGP` repository, the environment variables `FC` and
`FCFLAGS`, identifying the Fortran compiler and flags passed to it, need to be set. Here are some examples
used during development and testing.

To build any of the executables in `examples/` or `tests` the locations of the C and Fortran netCDF libraries
need to be set via environment variables `NCHOME` and `NFHOME`, and the variable `RRTMGP_ROOT` must be set to the
root of the RTE+RRTMGP installation.

## Gnu Fortran
`FC: gfortran-8` or `gfortran-9` or `gfortran-10`
### Debugging flags
`FCFLAGS: "-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -finit-real=nan -DUSE_CBOOL"`  
### Even stricter debugging flags
`FCFLAGS: "-ffree-line-length-none -m64 -std=f2008 -march=native -fbounds-check -fbacktrace -finit-real=nan -DUSE_CBOOL -pedantic -g -Wall"`  

## Intel Fortran
`FC: ifort`  
### Debugging flags
`FCFLAGS: "-m64 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132 -check bounds,uninit,pointers,stack -stand f08"`  
### Optimization flags:  
`FCFLAGS:"-m64 -O3 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132"`

## PGI Fortran
`FC: pgfortran` or `FC: nvfortran` (if using the Nvidia HPC SDK)
### Debugging flags
`FCFLAGS: "-g -Minfo -Mbounds -Mchkptr -Mstandard -Kieee -Mchkstk -Mallocatable=03  -Mpreprocess"`
### Optimization flags:  
`FCFLAGS: "-m64 -O3 -g -traceback -heap-arrays -assume realloc_lhs -extend-source 132"`
