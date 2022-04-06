---
project: RTE-Fortran
summary: RTE describes radiation problems in planetary atmospheres and computes radiative fluxes.
email: rrtmgp@aer.com
preprocessor: gfortran -E
display: public
sort: permission-alpha
graph: true
md_extensions: markdown.extensions.toc
title: RTE Fortran interfaces
src_dir: ../../rte
exclude_dir: ../../rte/kernels
exclude_dir: ../../rte/kernels-openacc
output_dir: ../../public/reference/rte-fortran-interface
...

These pages provide a programmer's view of the Fortran user interface to RTE.

RTE aspires to follow a set of coding conventions:

- Module names start with `mo_`, class/type names with `ty_`.
- Real variables are defined with working precision `wp`; logical variables with `wl`.
    Both Fortran KIND parameters are set in [mo_rte_kind.F90]
- Most procedures in RTE and RRTMGP are functions which return a string. A non-empty string indicates an error of some kind.
- RTE and RRTMGP operate on multiple columns (profiles) at once. Problems are dimensioned by column, layer,
    and spectral quadrature point.
- RTE and RRTMGP are agnostic to vertical ordering
- Units are MKS

The optical properties module provides a way to specify the spectrally-dependent
optical properties of the atmosphere consistent

- `init()` routines to specify e.g. the spectral discretization
- `load()` routines to provide data (e.g. lookup tables) needed to do a calculation
- `alloc()` routines to allocate memory once the problem size is known
- `finalize()` routines to reset variables to an un-initialized state
- Some classes have `get_ncol()` and `get_nlay()` methods to report the problem size
- Some classes have get_subset() methods to extract values along the column dimension


<!---
## How to Read This Documentation

Start with the [README] and the [Tutorial](./page/Tutorial.html).
Additionally, there is a page that provides a higher level organizational overview that you can find [here](./page/Organized_Listing.html).

The listings may not be exhaustive. Full listings are available using the links at the top of the page, where
there's also a search bar.

Take me back to the [User Documentation].

[README]: https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/README.md
[User Documentation]: ../../index.html
-->
[mo_rte_kind.F90]: ./module/mo_rte_kind.html
