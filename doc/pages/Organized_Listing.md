---
title: Organized Listing
---

[TOC]

## lists describing a providing links to specific items

The most involved class structure is in [mo_optical_properties.F90]. Here we have an  base class (`ty_optical_props`) that defines ~20 type-bound procedures, an abstract sub-class (`ty_optical_props_arry`) that defines three more deferred interfaces, and three sub-classes of ty_optical_props_arry (`ty_optical_props_1scl`, `ty_optical_props_2str`, `ty_optical_props_nstr` that add new procedures (e.g. `alloc_2str`)

[mo_optical_properties.F90]: ../sourcefile/mo_optical_props.F90.html