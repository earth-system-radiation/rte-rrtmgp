---
title: Organized Listing
---

[TOC]

## Describing and providing links to specific items

The most involved class structure is in [mo_optical_props.F90]. Here we have an  base class ([ty_optical_props]) that defines ~20 type-bound procedures, an abstract sub-class ([ty_optical_props_arry]) that defines three more deferred interfaces, and three sub-classes of [ty_optical_props_arry] ([ty_optical_props_1scl], [ty_optical_props_2str], [ty_optical_props_nstr] that add new procedures (e.g. [init_and_alloc_2str])

[mo_optical_props.F90]: ../module/mo_optical_props.html
[ty_optical_props]: ../type/ty_optical_props.html
[ty_optical_props_arry]: ../type/ty_optical_props_arry.html
[ty_optical_props_1scl]: ../type/ty_optical_props_1scl.html
[ty_optical_props_2str]: ../type/ty_optical_props_2str.html
[ty_optical_props_nstr]: ../type/ty_optical_props_nstr.html
[init_and_alloc_2str]: ../proc/init_and_alloc_2str.html