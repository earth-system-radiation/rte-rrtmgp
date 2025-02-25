---
layout: page
title: Working with optical property and source function values
---

# Shared behavior

Variables of class `ty_optical_props_arry` must be allocated as well as initialized before use. Because the argument list is slightly different among the sub-classes each allocation routine has a distinct name. Users must call one of `err_msg = op%alloc_1scl(ncol, nlay)`, `err_msg = op%alloc_2str(ncol, nlay)` or `err_msg = op%alloc_nstr(nmom, ncol, nlay)`.

Initialization from an existing variable of type `ty_optical_props` can be combined with allocation: `err_msg = op%alloc_2str(ncol, nlay, op1)`, etc.

Integer functions `op%get_ncol()` and `op%get_nlay()` return the number of columns or layers in the domain. Class `ty_optical_props_nstr` has a `get_nmom() `function.

`err_msg = op%subset(start, n, subset)` extracts columns `start` to `start+n-1` from variable `op` into variable `subset`. For optical properties assumptions are made if `op` and `subset` are of different classes.

Class `ty_source_func_lw` share these behaviors, though allocation and/or initialization does not require specifying the type (i.e. `err_msg = source_func%alloc(ncol, nlay)`).

# Working with optical property values

Calling `err_msg = op%delta_scale(for)` will delta-scale the optical properties to remove strong forward peaks in the phase function. Optional argument `for` is an real array of forward-scattering fractions with the same extents as `op%tau`. If this argument is omitted the assumption `for = g*g` is made.

Optical properties can be added together: `err_msg = op1%increment(op2)` adds the optical properties of op1 to those of op2 making reasonable assumptions if the two variables are not of the same class.

`err_msg = op%validate()` checks to be sure that all values are valid: optical depth not less than 0, single-scattering albedo between 0 and 1, asymmetry parameter between -1 and 1.

# Reference

A complete description of the [optical properties](rte-fortran-interface/type/ty_optical_props_arry.html) and [source function](./rte-fortran-interface/type/ty_source_func_lw.html) implementations is also generated from the source code.
