
#include "expand_and_transpose.h"

// Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
void expand_and_transpose(OpticalProps const &ops, real2d const &arr_in, real2d const &arr_out) {
  using yakl::intrinsics::size;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  int ncol  = yakl::intrinsics::size(arr_in, 2);
  int nband = ops.get_nband();
  int ngpt  = ops.get_ngpt();
  int2d limits = ops.get_band_lims_gpoint();
  // for (int iband=1; iband <= nband; iband++) {
  //   for (int icol=1; icol <= ncol; icol++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nband,ncol) , YAKL_LAMBDA (int iband, int icol) {
    for (int igpt=limits(1,iband); igpt <= limits(2,iband); igpt++) {
      arr_out(icol, igpt) = arr_in(iband,icol);
    }
  });
}

#ifdef RRTMGP_ENABLE_KOKKOS
// Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
void expand_and_transpose(OpticalPropsK const &ops, real2dk const &arr_in, real2dk const &arr_out)
{
  int ncol  = arr_in.extent(1);
  int nband = ops.get_nband();
  int ngpt  = ops.get_ngpt();
  int2dk limits = ops.get_band_lims_gpoint();
  // for (int iband=1; iband <= nband; iband++) {
  //   for (int icol=1; icol <= ncol; icol++) {
  Kokkos::parallel_for( MDRangeP<2>({0,0}, {nband,ncol}) , KOKKOS_LAMBDA (int iband, int icol) {
    for (int igpt=limits(0,iband); igpt <= limits(1,iband); igpt++) {
      arr_out(icol, igpt) = arr_in(iband,icol);
    }
  });

}
#endif
