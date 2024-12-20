
#pragma once

#include "rrtmgp_const.h"
#include "mo_optical_props.h"

#ifdef RRTMGP_ENABLE_YAKL
// Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
void expand_and_transpose(OpticalProps const &ops, real2d const &arr_in, real2d const &arr_out);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
// Expand from band to g-point dimension, transpose dimensions (nband, ncol) -> (ncol,ngpt)
template <typename RealT, typename LayoutT, typename DeviceT,
          typename ArrInT, typename ArrOutT>
void expand_and_transpose(OpticalPropsK<RealT, LayoutT, DeviceT> const &ops,
                          ArrInT const &arr_in, ArrOutT const &arr_out)
{
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  int ncol  = arr_in.extent(1);
  int nband = ops.get_nband();
  int ngpt  = ops.get_ngpt();
  Kokkos::View<int**, LayoutT, DeviceT> limits = ops.get_band_lims_gpoint();
  // for (int iband=1; iband <= nband; iband++) {
  //   for (int icol=1; icol <= ncol; icol++) {
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ncol, nband}) , KOKKOS_LAMBDA (int icol, int iband) {
    for (int igpt=limits(0,iband); igpt <= limits(1,iband); igpt++) {
      arr_out(icol, igpt) = arr_in(iband,icol);
    }
  }));
}
#endif
