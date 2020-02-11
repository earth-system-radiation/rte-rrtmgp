# All-sky example

This example provides a modestly more realistic setting the clear-sky problem done in the parallel `rfmip-clear-sky` directory in that it does an 'all-sky' calculation including both gases and clouds. It may be useful as a tutorial for users integrating RTE+RRTMGP into host models. We are also using it as a regression test for continuous integration and as a testbed for accelerators (currently GPUs using OpenACC).

The example uses the first of the Garand atmosphere used for developing RRTMGP, as described in the [paper](https://doi.org/10.1029/2019MS001621) documenting the code, repeats the column a user-specified number of times, computes the optical properties of an arbitrary cloud in each column, and computes the broadband fluxes. Fractional cloudiness is not considered, and the clouds are extensive but quite boring, with uniform condensate amount and particle size everywhere (though with different values for liquid and ice).

1. Build the RTE+RRTMGP libraries in `../../build/`. This will require setting environmental variables `FC` for the Fortran compiler and `FCFLAGS`, or creating `../../build/Makefile.conf` with that information.
2. Build the executables in this directory, which will require providing the locations of the netCDF C and Fortran libraries and module files as environmental variables (NCHOME and NFHOME) or via file `Makefile.libs`
4. Use Python script `run-rfmip-examples.py` to run the examples. The script takes some optional arguments, see `run-rfmip-examples.py -h`
5. Python script `compare-to-reference.py` will compare the results to reference answers for 128 columns, produced on a Mac with Intel 19 Fortran compiler.
