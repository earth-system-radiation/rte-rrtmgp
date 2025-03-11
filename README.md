[![Conda Version](https://img.shields.io/conda/vn/conda-forge/rte_rrtmgp.svg)](https://anaconda.org/conda-forge/rte_rrtmgp)
[![Conda Platforms](https://img.shields.io/conda/pn/conda-forge/rte_rrtmgp.svg)](https://anaconda.org/conda-forge/rte_rrtmgp)

The RTE+RRTMGP libraries and data files can be installed via `conda-forge`.

# Documentation

[RTE+RRTMGP's documentation](https://earth-system-radiation.github.io/rte-rrtmgp/) contains
a mix of automatically-generated pages and hand-written descriptions. The Wiki is deprecated.

# RTE+RRTMGP

This is the repository for RTE+RRTMGP, a set of codes for computing radiative fluxes in planetary atmospheres. RTE+RRTMGP is described in a [paper](https://doi.org/10.1029/2019MS001621) in [Journal of Advances in Modeling Earth Systems](http://james.agu.org).

RRTMGP uses a correlated _k_-distribution to provide an optical description (absorption and possibly Rayleigh optical depth) of the gaseous atmosphere, along with the relevant source functions, on a pre-determined spectral grid given temperatures, pressures, and gas concentration. The k-distribution currently distributed with this package is applicable to the Earth's atmosphere under present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source functions. The fluxes are normally summarized or reduced via a user extensible class.

## Building the libraries, examples, and unit-testing codes.

A description of building RTE+RRTMGP with `CMake` is described in the [documentation](https://earth-system-radiation.github.io/rte-rrtmgp/how-tos/).

See also the `autoconf` branch for a Gnu autotools build system.

## Examples

Two examples are provided in `examples/`, one for clear skies and one including clouds. Directory `tests/` contains regression testing (e.g. to ensure that answers are independent of orientation) and unit testing (to be sure all the code paths are tested). See the README file and codes in each directory for further information.

## Citing the code

Code releases are archived at Zenodo. All releases are available at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3403172.svg)](https://doi.org/10.5281/zenodo.3403172).
The current release is available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14841903.svg)](https://doi.org/10.5281/zenodo.14841903)

Please cite the code using these DOIs and the information in the `CITATION.cff` file in addition to the reference [paper](https://doi.org/10.1029/2019MS001621)

## Acknowledgements

The development of RTE+RRTMGP has been funded in the US by the Office of Naval Research, NASA, NOAA, and the Department of Energy. We
are grateful for contributions from a range of collaborators at institutions including the Swiss Supercomputing Center,
the German Climate Computing Center, and Nvidia.
