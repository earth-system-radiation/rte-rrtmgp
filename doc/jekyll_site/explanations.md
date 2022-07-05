---
layout: "page"
title: "Explanations"
---

# Explanations will live here

### Spectral discretization

The spectral properties of the atmosphere and the source functions depend on electromagnetic wavelength (or frequency or wavenumber). RTE treats this spectral dependence by dividing the spectrum into one or more _bands_, each of which represents a continuous set of wavelengths/frequencies. Bands may be further sub-divided into _g-points_ (the language is borrowed from _k_-distributions). Each _g_-point is treated as a independent psudo-monchromatic calculation but there is no inherent mapping between _g_-points and wavelengths; the sum over _g_-points is the band-average value.

The bands defined by RRTMGP cover the full spectrum of radiation emitted by the Sun and Earth: these are _broadband_ calculations. In RRTMGP the bands are continuous so that the ending wavelength of one band is the starting wavelength of the next.  
