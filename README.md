# Documentation

[RTE+RRTMGP's GitHub Pages site](https://earth-system-radiation.github.io/rte-rrtmgp/) contains
a mix of automatically-generated documentation and hand-written descriptions. The documentation is
incomplete and evolving. Thanks to the folks at [Sourcery Institute](https://www.sourceryinstitute.org)
for help in setting this up. 

For the moment the [Wiki](https://github.com/earth-system-radiation/rte-rrtmgp/wiki) may also be useful.

# RTE+RRTMGP

This is the repository for RTE+RRTMGP, a set of codes for computing radiative fluxes in planetary atmospheres. RTE+RRTMGP is described in a [paper](https://doi.org/10.1029/2019MS001621) in [Journal of Advances in Modeling Earth Systems](http://james.agu.org).

RRTMGP uses a k-distribution to provide an optical description (absorption and possibly Rayleigh optical depth) of the gaseous atmosphere, along with the relevant source functions, on a pre-determined spectral grid given temperatures, pressures, and gas concentration. The k-distribution currently distributed with this package is applicable to the Earth's atmosphere under present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source functions. The fluxes are normally summarized or reduced via a user extensible class.


## Building the libraries, examples, and unit-testing codes.

1. Set environment variables `FC` (the Fortran 2003 compiler) and `FCFLAGS` (compiler flags). Examples are provided in the `Compiler-flags.md` file.
2. Set environment variables `RRTMGP_ROOT` to the top-level RTE+RRTMGP directory and `RTE_KERNELS` to `openacc` if you want the OpenACC/OpenMP kernels rather than the default.
3. `make libs` in the top-level directory will make the RTE and RRTMGP libraries.
4. The examples and testing codes use netCDF. Set the variables `NCHOME` and `NFHOME` to the roots of the C and Fortran netCDF installations, then `make tests` to build and run these. (A few files need to be downloaded for `examples/rfmip-clear-sky`. The default is to download these with `wget` but a Python script is also available.)
5. Evaluating the results of the tests requires `Python` and the packages described in `environment.yml`. Comparisons can be made with `make check` in the top level directory.
6. `make` invoked without a target in the top level attempts all three steps.

## Examples

Two examples are provided in `examples/`, one for clear skies and one including clouds. Directory `tests/` contains regression testing (e.g. to ensure that answers are independent of orientation) and unit testing (to be sure all the code paths are tested). See the README file and codes in each directory for further information.
