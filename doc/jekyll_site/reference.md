---
layout: page
title: Reference (technical documentation)
---
# Reference/Technical documentation

RTE aspires to follow a set of coding conventions:

- Module names start with `mo_`, class/type names with `ty_`.
- Real variables are defined with working precision `wp`; logical variables with `wl`.
    Both Fortran KIND parameters are set in `mo_rte_kind.F90`
- Most procedures in RTE and RRTMGP are functions which return a string. A non-empty string indicates an error of some kind.
- RTE and RRTMGP operate on multiple columns (profiles) at once. Problems are dimensioned by column, layer,
    and spectral quadrature point.
- RTE and RRTMGP are agnostic to vertical ordering
- Units are MKS

## Fortran user-facing class interfaces

[RTE    Fortran interface](./reference/rte-fortran-interface/index.html)

[RRTMGP Fortran interface](./reference/rrtmgp-fortran-interface/index.html)

## Kernel interfaces

[RTE kernels:    Fortran interface](./reference/rte-kernels/index.html)

[RRTMGP kernels: Fortran interface](./reference/rrtmgp-kernels/index.html)
