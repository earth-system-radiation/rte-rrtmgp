---
layout: page
title: RTE and RRTMGP - optics
---

# General interfaces (procedure signatures) for gas optics

RTE includes the abstract class `ty_gas_optics` which defines the interfaces for computing the optical properties and radiation source functions for the gaseous atmosphere given pressure, temperature, and composition information. Both interfaces are accessed through `go%gas_optics()` (where `go` is a variable of class `ty_gas_optics`). Users can request output suitable for use in emission/absorption calculations or in scattering calculations (see the documentation of [optical properties](./rte-optical-props.html)). The Planck source function throughout the atmosphere is returned in a variable of `ty_source_func_lw` (as per [this description](./rte-optical-props.html)) while the top-of-atmosphere spectral flux is returned as an array. Requesting a source function for which the source function is unknown (i.e. asking for the top-of-atmosphere flux from a gas optics scheme for terrestrial radiation) returns an error.

# Describing compostion

Class `ty_gas_concs` provides a flexible Fortran-compatible way of describing gas concentrations. Gases are indexed by their chemical formula (i.e. "h2o") and may be provided as scalars, as 2D fields depending on column and level, or as 1D fields depending on level. Futher information is available as an [overview](../reference/gas-concentrations-overview.md).

# RRTMGP gas, cloud, and aerosol optics

Specific gas optics schemes such as [RRTMGP](../reference/rrtmgp-fortran-interface/index.html) implement the functions described in the interface. They may also implement other routines, for example to load data into the class at initialization. This allows scheme-specific code to be confined to initialization routines.

RRTMGP also implements schemes for computing the optical properties of [clouds](../reference/rrtmgp-fortran-interface/module/mo_cloud_optics_rrtmgp.html) and [aerosols](../reference/rrtmgp-fortran-interface/module/mo_aerosol_optics_rrtmgp_merra.html). Like the gas optics these schemes are initialized with data but the input variables are scheme-specific so there is no generic interface. See the [overview](../reference/rrtmgp-overview.html) for more details.
