---
project: RRTMGP-Fortran
summary: RRTMGP is a correlated k-distribution for computing fluxes in earth's atmosphere.
preprocessor: gfortran -E
display: public
sort: permission-alpha
graph: true
md_extensions: markdown.extensions.toc
author: The RTE+RRTTMGP consortium
github: https://github.com/earth-system-radiation/
license: by
title: RRTMGP Fortran interfaces
src_dir: ../../rrtmgp
exclude_dir: ../../rrtmgp/kernels
exclude_dir: ../../rrtmgp/kernels-openacc
output_dir: ../../public/reference/rrtmgp-fortran-interface
...
These pages provide a programmer's view of the Fortran user interface to RRTMGP.

RRTMGP provides a class [ty_gas_optics_rrtmgp](./type/ty_gas_optics_rrtmgp.html) that implements
the `gas_optics()` and other procedure(s) defined by the  [ty_gas_optics](./type/ty_gas_optics.html)
abstract class. The class is used to compute the spectrally-varying optical properties of the
gaseous atmosphere given temperature, pressure, and gas concentrations. Each instance of the
variable is ["loaded"](./type/ty_gas_optics_rrtmgp.html#boundprocedure-load) with data from netCDF
files in the `rrtmgp/data` directory. Depending on the data provided the variable can be used
or radiation emitted by the atmosphere and surface ("longwave") of for for radiation emitted
by the planet's star ("shortwave").

The class implements both the longwave/internal sources and
shortwave/external sources versions of the
[`gas_optics`](./type/ty_gas_optics_rrtmgp.html#boundprocedure-gas_optics~2) procedure.
The longwave version reports Planck sources at layer centers and layer interfaces (levels)
while the shortwave version reports the spectrally-varying stellar radiation
Calling the longwave routine (by providing the longwave-relevant arguments)
when the variable has been initialized with shortwave data triggers a run-time error.

The user interface uses the [ty_gas_concs](./module/mo_gas_concentrations.html) type
to represent the volume mixing ratios needed as input. Output suitable for
scattering emission, two-stream, or multi-stream calculations are provided
depending on which sub-class of RTE's
[ty_optical_props_arry](./rte-fortran-interface/module/mo_optical_props.html#type-ty_optical_props_arry)
are provided. Planck source functions, if requested, are reported in a variable
of type [ty_source_func_lw.](./rte-fortran-interface/type/ty_source_func_lw.html)

The listings below may not be exhaustive.
To see the full listings use the links at the top of the page.
There is a search bar in the top right.

Return to the [Documentation overview] or the [reference overview].

[Documentation overview]: ../../index.html
[reference overview]: ../../reference.html
