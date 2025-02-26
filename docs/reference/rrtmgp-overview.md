---
layout: page
title: Working with RRTMGP gas optics
---

**Initialization**:
Variables of class `ty_gas_optics_rrtmgp` require a set of lookup tables in order to compute optical properties from atmospheric state and composition. These data are provided in netCDF files in the [rrtmgp-data repository](https://github.com/earth-system-radiation/rrtmgp-data). The fields in the netCDF file must be passed to the `load()` functions, invoked as `err_msg = go%load(...)`, before the gas optics class can be used.

The longwave/shortwave variant is set depending on what data are provided to `load()`. Logical functions`go%source_is_external()` and `go%source_is_internal()` can be used to determine the variant.

RRTMGP does not read the netCDF files directly (we expect the code to be used in many environment with specific requirements) but a straightforward implementation is available in `examples/rfmip-clear-sky/`.

**Computing optimal angles for longwave radiative transfer**:
Longwave flux can be estimated from intensity calculations made at one or more angles; using more than one angle increases accuracy at the expense of more computation. RRTMGP gas optics can provide an empirical estimate, based on the atmospheric opacity, of the single spectrally-dependent angle at which to make such calculations to minimize error. For example:

```
    real(wp), dimension(ncol, ngpt) :: lw_Ds

    ! Compute the optical properties of the clear atmosphere in atmos
    call stop_on_err(k_dist%gas_optics(..., atmos)

    ! Compute the spectrally dependent optimal angle
    call stop_on_err(k_dist%compute_optimal_angles(atmos, lw_Ds))

    ! Use the optimal angle in flux calculations
    call stop_on_err(rte_lw(atmos, ..., lw_Ds=lw_Ds))
```

**Modifying the solar source function**:
The solar source function provided by the solar variant of `ty_gas_optics_rrtmgp` can be modified via calls. The total solar irradiance can be set directly with a call to `err_msg = go%set_tsi(tsi)`. Alternatively any or all of parameters of the NRLSSI2 model of solar variability (total solar irradiance, Bremen facular index, and sunspot index) can be provided via `err_msg = go%set_solar_variability(mg_index, sb_index, tsi)`. Values persist until they are explicitly reset. Directory `extensions/solar_variability/` contains code and data for determining these indices as a function of position in an average solar cycle.
