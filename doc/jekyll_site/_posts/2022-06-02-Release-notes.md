---
layout: post
title:  "2022-06-02: Release notes"
categories: Release-notes
---

Commit [a4fe30c](https://github.com/earth-system-radiation/rte-rrtmgp/commit/a4fe30cf4dab2e5fd8d3ab6f11683a82ae584475) 
to branch `main` makes the following changes:

- Solar zenith angle can vary with height, allowing for calculations on a pseudo-spherical earth
- New documentation site, partly hand-written (still essentially empty) and partly auto-generated from comments in code (partly complete).
- Workaround for Intel compilers on specific AMD processors
- Python dependencies made explicit in environment.yml files
- Move from wget to Python for downloading files in testing
- Public visibility and C-bindings restricted to just a few kernels
- Assorted consistency and prettifying
