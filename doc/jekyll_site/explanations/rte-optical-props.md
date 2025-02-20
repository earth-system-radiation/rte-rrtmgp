---
layout: page
title: RTE - optical properties
---

# Optical properties and spectral discretization

The spectral properties of the atmosphere and the source functions depend on electromagnetic wavelength (or frequency or wavenumber). RTE treats this spectral dependence by dividing the spectrum into one or more _bands_, each of which represents a continuous set of wavelengths/frequencies/wavenumbers. Bands may be further sub-divided into _g-points_ (the language is borrowed from _k_-distributions). Each _g_-point is treated as a independent psudo-monchromatic calculation but there is no inherent mapping between _g_-points and wavelengths; the sum over _g_-points is the band-average value.

Bands are defined by their bounds, expressed as wavenumbers in 1/cm, and by their staring and ending (integer) g-points. A spectral discretization defined only on bands is represented with one _g_-point per band. A set of monochromatic caclulations may be represented by additionally setting the upper and lower wavenumbers for each band to the same value.

Class `ty_optical_props` implements a range of [procedures](./reference/optical-props-overview.html) for representing and manipulating the spectral discretization. These capabilities are inherited by classes that represent arrays of optical properties (i.e. optical depth) or sources of radiation and by those used to compute optical properties from a physical description of the atmosphere.

# RRTMGP's spectral discretization

The bands defined by RRTMGP cover the full spectrum of radiation emitted by the Sun and Earth: these are _broadband_ calculations. In RRTMGP the bands are continuous so that the ending wavelength of one band is the starting wavelength of the next.
