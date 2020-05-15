
#include "expand_and_transpose.h"

// Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
void expand_and_transpose(OpticalProps const &ops, real2d const &arr_in, real2d &arr_out) {
  int ncol  = size(arr_in, 2);
  int nband = ops.get_nband();
  int ngpt  = ops.get_ngpt();
  int2d limits = ops.get_band_lims_gpoint();
  // for (int iband=1; iband <= nband; iband++) {
  //   for (int icol=1; icol <= ncol; icol++) {
  parallel_for( Bounds<2>(nband,ncol) , YAKL_LAMBDA (int iband, int icol) {
    for (int igpt=limits(1,iband); igpt <= limits(2,iband); igpt++) {
      arr_out(icol, igpt) = arr_in(iband,icol);
    }
  });
}


