This is the repository for RTE+RRTMGP, a set of codes for computing radiative fluxes in planetary atmospheres.

RRTMGP uses a k-distribution to provide an optical description (absorption and possibly Rayleigh optical depth) of the gaseous atmosphere, along with the relevant source functions, on a pre-determined spectral grid given temperatures, pressures, and gas concentration. The k-distribution currently distributed with this package is applicable to the Earth's atmosphere under present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source functions. The fluxes are normally summarized or reduced via a user extensible class.

Documentation and example programs to come relatively quickly. Meanwhile for questions please contact Robert Pincus and Eli Mlawer at rrtmgp@aer.com.

# Building the libraries.

1. `cd build`
2. Create a file `Makefile.conf` defining make variables `FC` (the Fortran 2003 compiler) and `FCFLAGS` (compiler flags). Alternately  link to an existing file or set these as environment variables.
3. `make`
