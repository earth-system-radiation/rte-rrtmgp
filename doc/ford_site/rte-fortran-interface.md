---
project: RTE-Fortran
summary: RTE describes radiation problems in planetary atmospheres and computes radiative fluxes.
preprocessor: gfortran -E
display: public
sort: permission-alpha
graph: true
md_extensions: markdown.extensions.toc
author: The RTE+RRTTMGP consortium
github: https://github.com/earth-system-radiation/
license: by
title: RTE Fortran interfaces
src_dir: ../../rte
exclude_dir: ../../rte/kernels
exclude_dir: ../../rte/kernels-openacc
output_dir: ../../public/reference/rte-fortran-interface
...

These pages provide a programmer's view of the Fortran user interface to RTE.

Procedures in the two solvers ([rte_lw](./proc/rte_lw.html), [rte_sw](./proc/rte_sw.html))
require problems specified as sets of [optical properties](./module/mo_optical_props.html)
along with boundary conditions and/or [internal sources](./module/mo_source_functions.html) of radiation.
The reduction of spectrally- and spatially-detailed calculations that specifies the output
is defined in a [flux output type](./module/mo_fluxes.html)

The optical properties module provides a way to specify the spectrally-dependent
optical properties of the atmosphere consistent with calculations neglecting scattering,
using two-stream scattering, or using n-stream scattering.

- `init()` routines to specify e.g. the spectral discretization
- `load()` routines to provide data (e.g. lookup tables) needed to do a calculation
- `alloc()` routines to allocate memory once the problem size is known
- `finalize()` routines to reset variables to an un-initialized state
- Some classes have `get_ncol()` and `get_nlay()` methods to report the problem size
- Some classes have `get_subset()` methods to extract values along the column dimension

The listings below may not be exhaustive.
To see the full listings use the links at the top of the page.
There is a search bar in the top right.

Return to the [Documentation overview] or the [reference overview].

[Documentation overview]: ../../index.html
[reference overview]: ../../reference.html

<!---
## How to Read This Documentation

Start with the [README] and the [Tutorial](./page/Tutorial.html).
Additionally, there is a page that provides a higher level organizational overview that you can find [here](./page/Organized_Listing.html).

The listings may not be exhaustive. Full listings are available using the links at the top of the page, where
there's also a search bar.

Take me back to the [User Documentation].

[README]: https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/README.md
[User Documentation]: ../../index.html
[mo_rte_kind.F90]: ./module/mo_rte_kind.html
-->
