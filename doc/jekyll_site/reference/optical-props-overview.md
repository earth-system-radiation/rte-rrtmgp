---
layout: page
title: Working with optical properties and spectral discretizations
---

This page provides an overview of methods available for working with spectral discretizations as used by the optical properties class and its descendents. Further details are available in the [auto-generated documentation](./reference/rte-fortran-interface/type/ty_optical_props.html).

Optical properties are defined by their spectral dependence which is described during [initialization](./reference/rte-fortran-interface/type/ty_optical_props.html).

The spectral dependence can be described via arrays: array `band_lims_wvn` with extents (2, number-of-bands) describes the number of bands beginning and ending wavenumber of each band (in MKS units i.e. inverse meters). Optional argument `band_lims_gpt` describes the beginning and ending g-point of each band; if this array isn't provided it's assumed values are available by band. An optional `name` may be useful in debugging.

Alternatively variables of type `ty_optical_props` can copy the spectral discretization from a previously-initialized variable e.g. `err_msg = op2%init(op1)`.

Real-valued functions ``` op%get_band_lims_wavenumber()`` and  ```op%get_band_lims_wavelength()\` return arrays of extent (2, number-of-bands) containing the spectral limits of each band expressed in wavelength or wavenumber = 1/wavelength (again in MKS units).

Procedure `op%is_initialized()` returns a logical value indicating if the spectral discretization has been provided. `op%get_nband()` and `op%get_ngpt()` return the (integer) number of band and g-points; if these numbers are the same the optical properties are defined by band.

Integer functions `op%convert_band2gpt()` and `op%convert_gpt2band()` provide a map between individual bands and g-points. Function `op%get_gpoint_bands()` returns an integer of length op%get_ngpt() with the band number to which each g-point belongs (i.e. it's the vector version of `op%convert_gpt2band()`).

Comparison functions `op1%bands_are_equal(op2)` and `op1%gpoints_are_equal(op2)` return logical values indicating whether the two sets of optical properties share the same band or g-point discretization.
