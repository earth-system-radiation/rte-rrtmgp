# rfmip-rrtmgp
This directory contains programs and support infrastructure for running
the [RTE+RRTMGP](https://github.com/RobertPincus/rte-rrtmgp) radiation parameterization for the
[RFMIP](https://www.earthsystemcog.org/projects/rfmip/) cases.

1. Build the RTE+RRTMGP libraries in `../../build/`. This will require setting
environmental variables `FC` for the Fortran compiler and `FCFLAGS`, or creating
`../../build/Makefile.conf` with that information.
2. Build the executables in this directory, which will require providing the
locations of the netCDF C and Fortran libraries and module files as environmental
variables (NCHOME and NFHOME) or via file `Makefile.libs`
3. Use Python script `stage_files.py` to download relevant files from the
[RFMIP web site](https://www.earthsystemcog.org/projects/rfmip/resources/).This script invokes another Python script to create empty output files.
4. Use Python script `run-rfmip-examples.py` to run the examples. The script takes
some optional arguments, see `run-rfmip-examples.py -h`
5. Python script `compare-to-reference.py` will compare the results to reference
answers produced on a Mac with Intel 19 Fortran compiler. Differences are normally
within 10<sup>-6</sup> W/m<sup>2</sup>.

The Python scripts require modules `netCDF4`, `numpy`, and `xarray`.
