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

A description of building RTE+RRTMGP with an ad hoc homemade system is described in the [documentation](https://earth-system-radiation.github.io/rte-rrtmgp/how-tos/).

See also the `autoconf` branch for a Gnu autotools build system. 

## Examples

Two examples are provided in `examples/`, one for clear skies and one including clouds. Directory `tests/` contains regression testing (e.g. to ensure that answers are independent of orientation) and unit testing (to be sure all the code paths are tested). See the README file and codes in each directory for further information.

## Citing the code 

Code releases are archived at Zenodo. All releases are available at 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403172.svg)](https://doi.org/10.5281/zenodo.3403172). 
The current release is available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7521518.svg)](https://doi.org/10.5281/zenodo.10211873)

Please cite the code using these DOIs and the information in the `CITATION.cff` file in addition to the reference [paper](https://doi.org/10.1029/2019MS001621)

## Acknowledgements

The development of RTE+RRTMGP has been funded in the US by the Office of Naval Research, NASA, NOAA, and the Department of Energy. We 
are grateful for contributions from a range of collaborators at institutions including the Swiss Supercomputing Center, 
the German Climate Computing Center, and Nvidia. 