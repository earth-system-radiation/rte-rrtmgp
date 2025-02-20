---
layout: page
title: RTE - fluxes
---

RTE solves the radiative transfer equation for each spectral point independently but this detailed information typically isn't useful. Class `ty_fluxes` in module [mo_fluxes](./reference/rte-fortran-interface/module/mo_fluxes.html) in the RTE library provides a way to reduce the highly-detailed information based on precisely what the user needs. Class `ty_fluxes_broadband` provides an example implementation that reports broadband flux i.e. the sum over all spectral points at all levels.

# Fluxes interface

`[ty_fluxes](https://earth-system-radiation.github.io/rte-rrtmgp/reference/rte-fortran-interface/type/ty_fluxes.html)` is an abstract class that defines [two interfaces](./reference/rte-fortran-interface/type/ty_fluxes.html) defining type-bound functions. Function `reduce` is called within the RTE solvers. The input arguments include the spectrally-resolved fluxes up and down (and, optically, the direct-beam flux) and the spectral discretization. Logical function `are_desired` is called by RTE to check if the results of a calculation will be used.

Class `ty_fluxes` is abstract; it defines only the interfaces to these routines. Implementation is deferred to user classes that extend this class.

# Broadband fluxes

Class [ty_fluxes_broadband](./reference/rte-fortran-interface/type/ty_fluxes_broadband.html) in the same [module](./reference/rte-fortran-interface/module/mo_fluxes.html) provides an example of how to extend `ty_fluxes`. This class reports broadband values i.e. the sum over every element in the spectral dimensions. The class contains [four data fields](./reference/rte-fortran-interface/type/ty_fluxes_broadband.html). Before calling one of the RTE solvers with this class as output, users either allocate one or more of these output fields or point to an existing array. The fields must have the same number of columns and layers as the problem being solved. In class `ty_fluxes_broadband` logical function \`are_desired()\`\` returns true if one or more of the data fields points to storage and false otherwise.

# Other extensions

Class `ty_fluxes_byband` in subdirectory `extensions/` provides another example that makes use of the spectral discretization information.

It's useful to extend the class when particular subsets of fluxes, such as photosynthetically-active radiation at the surface, are required.
