---
layout: page
title: RTE - solvers
---

The [rte_lw()](./reference/rte-fortran-interface/proc/rte_lw.html) and [rte_sw()](./reference/rte-fortran-interface/proc/rte_sw.html) routines in the RTE library compute radiative transfer on a user-specified problem. Both routines allow for an upper boundary condition for spectrally-resolved diffuse radiation. `rte_sw()` is overloaded to allow users to specify either a single value of `mu0` per column or a value at each layer of each column.

If optical properties are specified as variables of type ty_optical_props_1scl (meaning extinction optical depth alone) in the longwave, radiative transfer is computed accounting for only emission and absorption. If optical properties are specified as variables of type `ty_optical_props_2str` (optical depth, single-scattering albedo, and asymmetry parameter) the rescaling and refinement method of [Tang et al. 2018](https://doi.org/10.1175/JAS-D-18-0014.1) is used by default; two-stream and adding methods maybe be chosen by setting optional argument `use_2stream` to `.TRUE.`

The sensitivity of broadband flux (W/m$^2$-K) to changes in surface temperature is available in optional output argument `flux_up_Jac` . Jacobians are not available if `use_2stream` is \`\`.TRUE.\`
