
#pragma once

#include "rrtmgp_const.h"
#include "mo_rrtmgp_util_reorder_kernels.h"

#ifdef RRTMGP_ENABLE_YAKL
// (x,y,z) -> (z,x,y)
void reorder123x312(int d1, int d2, int d3, real3d const &array, real3d const &array_out);

// (x,y,z) -> (z,y,x)
void reorder123x321(int d1, int d2, int d3, real3d const &array, real3d const &array_out);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
// (x,y,z) -> (z,x,y)
template <typename ArrayT, typename ArrayOutT>
void reorder123x312(int d1, int d2, int d3, ArrayT const &array, ArrayOutT const &array_out) {
  reorder_123x312_kernel(d1, d2, d3, array, array_out);
  std::cout << "WARNING: THIS IS NOT TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}

// (x,y,z) -> (z,y,x)
template <typename ArrayT, typename ArrayOutT>
void reorder123x321(int d1, int d2, int d3, ArrayT const &array, ArrayOutT const &array_out) {
  reorder_123x321_kernel(d1, d2, d3, array, array_out);
}
#endif
