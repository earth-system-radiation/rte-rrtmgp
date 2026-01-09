---
project: SSM optics
summary: RTE describes radiation problems in planetary atmospheres and computes radiative fluxes.
preprocess: false
display: public
sort: permission-alpha
graph: true
md_extensions: markdown.extensions.toc
author: The RTE+RRTTMGP consortium
github: https://github.com/earth-system-radiation/
license: by
output_dir: ../reference/ssm
src_dir: ../../ssm
---

These pages document the Simple Spectral Models of gas and cloud optics used by RTE-SSM.

The gas optics is described in [Williams (2026)](https://arxiv.org/abs/2508.09353) and [PR #379](https://github.com/earth-system-radiation/rte-rrtmgp/pull/379), and includes simple analytic representations to the absorption coefficients of H$_{2}$O and CO$_{2}$ line absorption (in the longwave) and H$_{2}$O and O$_{3}$ line absorption (in the shortwave).

The cloud optics is highly simplified and just prescribes a single mass absorption coefficient ($\kappa$), single scattering albedo ($\omega$), and asymmetry factor ($g$) for clouds in the longwave and shortwave, respectively.

The listings below are not exhaustive.
To see the full listings use the links at the top of the page.
There is a search bar in the top right.

Return to the [Documentation overview] or the [reference overview].

[documentation overview]: ../../index.html
[reference overview]: ../index.html
