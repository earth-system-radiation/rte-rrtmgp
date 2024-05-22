
#pragma once

#include "rrtmgp_const.h"
#include "mo_optical_props.h"

#ifdef RRTMGP_ENABLE_YAKL
// Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
void expand_and_transpose(OpticalProps const &ops, real2d const &arr_in, real2d const &arr_out);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
// Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
void expand_and_transpose(OpticalPropsK<> const &ops, real2dk const &arr_in, real2dk const &arr_out);
#endif
