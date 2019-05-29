# RTE+RRTMGP

This is the repository for RTE+RRTMGP, a set of codes for computing radiative fluxes in planetary atmospheres. RTE+RRTMGP is described in a [manuscript revised May 10, 2019](https://doi.org/10.1002/essoar.10500964.1) to [Journal of Advances in Modeling Earth Systems](http://james.agu.org). 

RRTMGP uses a k-distribution to provide an optical description (absorption and possibly Rayleigh optical depth) of the gaseous atmosphere, along with the relevant source functions, on a pre-determined spectral grid given temperatures, pressures, and gas concentration. The k-distribution currently distributed with this package is applicable to the Earth's atmosphere under present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source functions. The fluxes are normally summarized or reduced via a user extensible class.

Example programs and documenation are evolving - please see examples/ in the repo and Wiki on the project's Github page. Suggestions are welcome. Meanwhile for questions please contact Robert Pincus and Eli Mlawer at rrtmgp@aer.com.

## Building the libraries.

1. `cd build`
2. Set environment variables `FC` (the Fortran 2003 compiler) and `FCFLAGS` (compiler flags). Alternately create a Makefile.conf that sets these variables. You could also link to an existing file. 
3. Set environment variable `RTE_KERNELS` to `openacc` if you want the OpenACC kernels rather than the default. 
4. `make`

## Building and running the examples.

1. From the root RTE+RRTMGP directory: `cd examples/rfmip-clear-sky`. 
2. Set environment variables `NCHOME` and `NFHOME` to the root directories of the Netcdf C and Fortran libraries respectively. Set environment variable `RRTMGP_DIR` to the location of the libraries (`../../build`) in the default layout). 
3. `make`
4. Python scripts are provided to stage the files needed (`stage_files.py`),  run the examples `run-rfmip-examples.py`), and compare to results computed on an example host (`compare-to-reference.py`). The python scripts require modules xarray and netCDF.  
