---
project: rte-rrtmgp
summary: A set of codes for computing radiative fluxes in planetary atmospheres.
project_website: https://github.com/earth-system-radiation/rte-rrtmgp.git
email: rrtmgp@aer.com
src_dir: ../rrtmgp
src_dir: ../rte
src_dir: ../extensions
exclude: mo_gas_optics_kernels.F90
exclude: mo_rrtmgp_util_reorder_kernels.F90
exclude: mo_gas_optics_kernels.F90
exclude: mo_fluxes_broadband_kernels.F90
exclude: mo_optical_props_kernels.F90
exclude: mo_rte_solver_kernels.F90
exclude: mo_optical_props_kernels.F90
exclude: mo_rte_solver_kernels.F90
page_dir: ../doc/pages
preprocessor: gfortran -E
display: public
sort: permission-alpha
output_dir: ../public/reference/extensions/
graph: true
md_extensions: markdown.extensions.toc
...

Welcome to the rte-rrtmgp developer documentation.

## How to Read This Documentation

Start with the [README] and the [Tutorial](./page/Tutorial.html).
Additionally, there is a page that provides a higher level organizational overview that you can find [here](./page/Organized_Listing.html).

The listings below are not exhaustive.
To see the full listings use the links at the top of the page.
Also, if you know what you're looking for, there is a search bar in the top right.

[README]: https://github.com/earth-system-radiation/rte-rrtmgp/blob/main/README.md
