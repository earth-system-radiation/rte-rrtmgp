#include "mo_fluxes_byband_kernels.h"

#ifdef RRTMGP_ENABLE_YAKL
// Spectral reduction over all points
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &bnd_lims, real3d const &spectral_flux, real3d &byband_flux) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlev,ncol) , YAKL_LAMBDA (int ibnd, int ilev, int icol) {
    real bb_flux_s = 0.0;
    for (int igpt=bnd_lims(1,ibnd); igpt<=bnd_lims(2,ibnd); igpt++) {
      bb_flux_s += spectral_flux(icol,ilev,igpt);
    }
    byband_flux(icol,ilev,ibnd) = bb_flux_s;
  }));
}
// Compute net flux
void net_byband(int ncol, int nlev, int nbnd, real3d const &bnd_flux_dn, real3d const &bnd_flux_up, real3d &bnd_flux_net) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlev,ncol), YAKL_LAMBDA(int ibnd, int ilev, int icol) {
      bnd_flux_net(icol,ilev,ibnd) = bnd_flux_dn(icol,ilev,ibnd) - bnd_flux_up(icol,ilev,ibnd);
  }));
}
#endif

