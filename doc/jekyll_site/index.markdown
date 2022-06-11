---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults
layout: "home"
title: "RTE+RRTMGP documentation"
---
This is the documentation for RTE+RRTMGP, a set of codes for computing radiative
fluxes in planetary atmospheres. RTE+RRTMGP is described in a
[paper](https://doi.org/10.1029/2019MS001621) in
[Journal of Advances in Modeling Earth Systems](http://james.agu.org).

RRTMGP uses a k-distribution to provide an optical description (absorption and
  possibly Rayleigh optical depth) of the gaseous atmosphere, along with the
  relevant source functions, on a pre-determined spectral grid given
  temperatures, pressures, and gas concentration. The k-distribution currently
  distributed with this package is applicable to the Earth's atmosphere under
  present-day, pre-industrial, and 4xCO2 conditions.

RTE computes fluxes given spectrally-resolved optical descriptions and source
functions. The fluxes are normally summarized or reduced via a user extensible class.

# Work in progress

Comprehensive documentation for RTE and RRTMGP is very much under development.

We are planning to follow the [Diataxis](https://diataxis.fr/) framework
and produce [tutorials](tutorials/index.md), [how-to guides](howtos/index.html),
[explanations](explanations/index.html), and [reference](reference/index.html).

We are starting with the [reference documentation](reference/index.html),
auto-generated from the code itself. This is provided separately for
RTE and RRTMGP and for the user-facing classes and underlying computational kernels.
