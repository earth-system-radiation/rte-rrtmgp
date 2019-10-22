# RTE+RRTMGP

This is the repository for RTE+RRTMGP, a set of codes for computing radiative fluxes in planetary atmospheres. RTE+RRTMGP is described in a [paper](https://doi.org/10.1029/2019MS001621) in [Journal of Advances in Modeling Earth Systems](http://james.agu.org).

RRTMGP uses a k-distribution to provide an optical description (absorption and possibly Rayleigh optical depth) of the gaseous atmosphere, along with the relevant source functions, on a pre-determined spectral grid given temperatures, pressures, and gas concentration. The k-distribution currently distributed with this package is applicable to the Earth's atmosphere under present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source functions. The fluxes are normally summarized or reduced via a user extensible class.

Example programs and documenation are evolving - please see examples/ in the repo and Wiki on the project's Github page. Suggestions are welcome. Meanwhile for questions please contact Robert Pincus and Eli Mlawer at rrtmgp@aer.com.

## Building the libraries.

1. `cd build`
2. Set environment variables `FC` (the Fortran 2003 compiler) and `FCFLAGS` (compiler flags). Alternately create a Makefile.conf that sets these variables. You could also link to an existing file.
3. Set environment variable `RTE_KERNELS` to `openacc` if you want the OpenACC kernels rather than the default.
4. `make`

## Examples

Two examples are provided, one for clear skies and one including clouds. See the README file and codes in each directory for further information. 
