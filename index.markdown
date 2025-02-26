---
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults
layout: home
title: RTE+RRTMGP documentation
---

This is the documentation for RTE+RRTMGP, a set of codes for computing radiative
fluxes in planetary atmospheres. The initial release of RTE+RRTMGP is described in a
[paper](https://doi.org/10.1029/2019MS001621).
The code itself can be cited as
doi:[10.5281/zenodo.3403172](https://doi.org/10.5281/zenodo.3403172) or via the
DOI attached to each release.

RRTMGP uses a k-distribution to provide an optical description (absorption and
possibly Rayleigh optical depth) of the gaseous atmosphere, along with the
relevant source functions, on a pre-determined spectral grid given
temperatures, pressures, and gas concentration. The k-distribution currently
distributed with this package is applicable to the Earth's atmosphere under
present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source
functions. The fluxes are normally summarized or reduced via a user extensible class.

RTE+RRTMGP is a collaboration between Columbia University and Atmospheric and Environmental Research, Inc. with support from a growing community including the German Climate Computing Center, Nvidia, and the Swiss Supercomputing Center. Please open an [issue](https://github.com/earth-system-radiation/rte-rrtmgp/issues) with questions.

# Documentation is always a work in progress

We follow the [Diataxis](https://diataxis.fr/) framework
and produce [tutorials](./tutorials/index.html), [how-to guides](./how-tos/index.html),
[explanations](./explanations/index.html), and [reference](./reference/index.html).

Much of the [reference documentation](./reference/index.html) is
auto-generated from the code itself. This is provided separately for
RTE and RRTMGP and for the user-facing classes and underlying computational kernels.

We welcome contributions to the documentation via pull requests against the `main` branch
of the Github repository.
