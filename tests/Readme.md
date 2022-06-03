# Tests

This directory contains testing codes used in development.

## Verification and validation

The Makefile target `clear_sky_regression` is used for verification, meaning
ensuring the package produces results as expected (e.g. that the
rescaling approximation used for LW scattering reduces to the no-scattering
solution in the absence of clouds), and for validation, meaning 
comparisons among different approaches (e.g. different numbers of streams).

## Feature testing

The Makefile target `test_zenith_angle_spherical_correction` tests features related
to computing solar radiation in spherical geometry, including codes for computing
solar zenith angle as a function of height (`extensions/mo_zenith_angle_spherical_correction.F90`)
and fluxes computed with height-dependent zenith angle.

Module `mo_rcemip_profiles.F90` provides profiles of temperature, height, and gas
concentrations given an arbitrary set of pressure levels, based on the
[protocol](https://doi.org/10.5194/gmd-11-793-2018) for the
[Radiative-Convective Equilibrium MIP](https://myweb.fsu.edu/awing/rcemip.html).
