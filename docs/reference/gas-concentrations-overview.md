---
layout: page
title: Working with gas concentrations
---

Class ty_gas_concs from module mo_gas_concentrations in the RRTMGP library is the means for specifying gas concentrations to ty_gas_optics%gas_optics. The class has four user-facing procedures. The following examples assume the declaration `type(ty_gas_concs) :: gas_concs`

`err_msg = gas_optics%init(gas_names)` takes a vector of strings specifying the names of the gases that will be provided. Gas names are normally the chemical formula. Some compilers expect each string to be the same length, so shorter names may have to be padded with trailing spaces.

`err_msg = gas_optics%init(gas_name, value)`specifies the absolute volume mixing ratio for the named gas. Gas names are normally the chemical formula. Values may be specified as scalars, one-dimensional profiles assumed to apply layer-by-layer to all columns, or as a 2D field dimensioned (columns, layers). The number of columns is set the first time a 1D or 2D field is supplied; the number columns is determined from the first 2D field. Specifying a concentration for a gas not provided during initialization is an error.

`err_msg = gas_concs%get_subset(2, n, gas_conc_subset)` extracts columns start to start+n-1 from the class. Example:

``` call gas_concs%reset()`` will return variable  ```gas_concs\` to an un-initialized state. This is useful if, for example, the problem size has changed since concentrations were last specified.

`gas_concs%get_vmr()` is used within `ty_gas_optics%gas_optics()` to extract the gas concentrations. Asking for the vmr of a gas that wasn't provided during `init()` or hasn't been set through a call to `set_vmr()` is an error.
