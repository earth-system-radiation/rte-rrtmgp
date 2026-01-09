---
layout: page
title: Conventions, technical documentation
---

# Overviews

- [Optical properties and spectral discretization](./optical-props-overview.html)
- [Problem descriptions](./optical-props-src-funcs.html)
- [Gas concentrations](./gas-concentrations-overview.html)
- [Gas optics: generic interface, RRTMGP implementation](./rte-optics.md)

# Reading the code

RTE and its related optics packages aspire to follow a set of coding conventions:

- Module names start with `mo_`, class/type names with `ty_`.
- Real variables are defined with working precision `wp`; logical variables with `wl`.
  Both Fortran KIND parameters are set in `mo_rte_kind.F90`
- Most procedures are functions which return a string. A non-empty string indicates an error of some kind.
- Procedures operate on multiple columns (profiles) at once. Problems are dimensioned by column, layer,
  and spectral quadrature point.
- Procedures are agnostic to vertical ordering
- Units are MKS
- Procedures (with the exception of testing code) do not perform I/O

# Radiative tranfer for energetics

[RTE Fortran interface](./rte-fortran-interface/index.html)

[RTE kernels: Fortran interface](./rte-kernels/index.html)

# RRTMGP optics

[RRTMGP Fortran interface](./rrtmgp-fortran-interface/index.html)

[RRTMGP kernels: Fortran interface](./rrtmgp-kernels/index.html)

# Simple spectral model 

[SSM Fortran interface and kernels](./ssm/index.html)
