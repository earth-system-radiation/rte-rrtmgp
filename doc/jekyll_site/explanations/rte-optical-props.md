---
layout: page
title: RTE - optical properties
---

# Optical properties and spectral discretization

The spectral properties of the atmosphere and the source functions depend on electromagnetic wavelength (or frequency or wavenumber). RTE treats this spectral dependence by dividing the spectrum into one or more _bands_, each of which represents a continuous set of wavelengths/frequencies/wavenumbers. Bands may be further sub-divided into _g-points_ (the language is borrowed from _k_-distributions). Each _g_-point is treated as a independent psudo-monchromatic calculation but there is no inherent mapping between _g_-points and wavelengths; the sum over _g_-points is the band-average value.

Bands are defined by their bounds, expressed as wavenumbers in 1/cm, and by their staring and ending (integer) g-points. A spectral discretization defined only on bands is represented with one _g_-point per band. A set of monochromatic caclulations may be represented by additionally setting the upper and lower wavenumbers for each band to the same value.

Class `ty_optical_props` implements a range of [procedures](../reference/optical-props-overview.html) for representing and manipulating the spectral discretization. These capabilities are inherited by classes that represent arrays of optical properties (i.e. optical depth) or sources of radiation and by those used to compute optical properties from a physical description of the atmosphere.

# Concrete sets of optical properties

RTE supports three concrete descriptions of radiative transfer problems expressed as arrays of discrete values on a domain with dimensions column, layer, spectral location (g-point).

- For problems considering only emission and absorption class `ty_optical_props_1scl` contains the data field `tau` representing absorption optical depth.
- For low-order problems considering scattering class `ty_optical_props_2str` contains the data fields `tau` (extinction optical depth), `ssa` (single-scattering albedo, and `g` (asymmetry parameter).
- For problems requiring higher-order treatments of scattering class `ty_optical_props_nstr` contains the data fields `tau`, `ssa`, and `p`. The latter has leading dimension `nmom` to accommodate the Legendre moments of a phase function. RTE does not contain a solution method for this class.

Each of these descriptions is a sub-class of the generic class `ty_optical_props_arry`. These classes share a range of [behaviors](../reference/optical-props-src-funcs.html) in addition to the procedures for working with spectral discretizations.

# (Longwave) source functions

Class `ty_source_func_lw` describes longwave source functions. It shares some behavior with `ty_optical_props_arry` but the data fields are not user-accessible.

# Digression: RRTMGP's spectral discretization

The bands defined by RRTMGP cover the full spectrum of radiation emitted by the Sun and Earth: these are _broadband_ calculations. In RRTMGP the bands are continuous so that the ending wavelength of one band is the starting wavelength of the next.
