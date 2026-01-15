---
layout: page
title: Building and testing
---

How to build the libraries, tests, and examples, run the tests, and verify the results

# In a nutshell

RTE+RRTMGP uses `CMake`. To generate a configuration using the `ninja` build system:
`cmake -S . -B build -G "Ninja"  -DCMAKE_BUILD_TYPE=RelWithDebInfo -DRTE_BUILD_TESTING=ON`

To build the code once a configuation has been generated:
`cmake --build build`

To run the examples and tests and check the results:
`ctest --test-dir build`

Evaluating the results of the tests requires `Python` and the packages described in `environment*.yml`.

See also possible values of [compiler flags](../reference/compiler-flags.html).

# Building and testing using (Gnu) autotools

Sergey Kosukhin and his colleagues at the Max Planck Institute for Meteorology
maintain the `autoconf` branch which adds Gnu `autotools` building to `main` branch.

# Supplying data

Data files needed for RRTMGP are available via a [data repository](https://github.com/earth-system-radiation/rrtmgp-data) or
as a [Zenodo archive](https://doi.org/10.5281/zenodo.7988260).

`ctest` fetches a specific version of this data for running the tests and examples.
