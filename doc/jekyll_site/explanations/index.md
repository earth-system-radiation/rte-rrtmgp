---
layout: page
title: Explanations
---

RTE+RRTMGP is a set of codes for computing radiative fluxes in planetary atmospheres. These pages provide a orientation to the code's architcture.

RTE [computes](./explanations/rte-solvers.html) radiative [_fluxes_](./explanations/rte-fluxes.html) given values of [_optical properties_](./explanations/rte-optical-props.html) , source functions, and boundary conditions.

RRTMGP computes the optical properties and source functions of the gaseous atmosphere given the distribution of temperature, pressure, and gas concentrations within the atmosphere.

Each of the italicized phrases above corresponds to a class in the Fortran 2003 interface which bundles some combination of code and data.
