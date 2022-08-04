---
layout: "page"
title: "Tutorials"
---

# Tutorials will live here

### Clear-sky calculation
A typical workflow for a clear-sky calculation is to
1. allocate memory
2. set gas concentration values
3. compute the optical properties of the gaseous atmosphere
4. compute radiative fluxes   

This repository contains all the pieces needed to perform a clear-sky calculation. An [example](https://github.com/RobertPincus/rte-rrtmgp/tree/master/examples/rfmip-clear-sky) is provided.

### All-sky calculation
An all-sky calculation is a small variant
1. allocate memory
2. set gas concentration values
3. compute the optical properties of the gaseous atmosphere
4. compute the optical properties of aerosol and add these to the optical properties
5. compute the optical properties of clouds and add these to the optical properties
6. compute radiative fluxes

An [example](https://github.com/earth-system-radiation/rte-rrtmgp/tree/main/examples/all-sky) of this workflow is available in the repository. The example also demonstrates how to complete an end-to-end calculation on a GPU using OpenACC. Users must provide methods for computing the optical properties of clouds and gases (an [example cloud optics class](https://github.com/earth-system-radiation/rte-rrtmgp/tree/main/extensions/cloud_optics)  is provided).
