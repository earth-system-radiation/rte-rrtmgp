
#include "mo_rrtmgp_util_reorder.h"

// (x,y,z) -> (z,x,y)
void reorder123x312(int d1, int d2, int d3, real3d const &array, real3d const &array_out) {
  reorder_123x312_kernel(d1, d2, d3, array, array_out);
  std::cout << "WARNING: THIS IS NOT TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// (x,y,z) -> (z,y,x)
void reorder123x321(int d1, int d2, int d3, real3d const &array, real3d const &array_out) {
  reorder_123x321_kernel(d1, d2, d3, array, array_out);
}

#ifdef RRTMGP_ENABLE_KOKKOS
void reorder123x312(int d1, int d2, int d3, real3dk const &array, real3dk const &array_out) {
  reorder_123x312_kernel(d1, d2, d3, array, array_out);
  std::cout << "WARNING: THIS IS NOT TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}



// (x,y,z) -> (z,y,x)
void reorder123x321(int d1, int d2, int d3, real3dk const &array, real3dk const &array_out) {
  reorder_123x321_kernel(d1, d2, d3, array, array_out);
}
#endif
