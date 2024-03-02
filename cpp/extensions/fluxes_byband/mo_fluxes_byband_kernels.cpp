#include "mo_fluxes_byband_kernels.h"

// Spectral reduction over all points
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &bnd_lims, real3d const &spectral_flux, real3d &byband_flux) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlev,ncol) , YAKL_LAMBDA (int ibnd, int ilev, int icol) {
    real bb_flux_s = 0.0;
    for (int igpt=bnd_lims(1,ibnd); igpt<=bnd_lims(2,ibnd); igpt++) {
      bb_flux_s += spectral_flux(icol,ilev,igpt);
    }
    byband_flux(icol,ilev,ibnd) = bb_flux_s;
  });
}
// Compute net flux
void net_byband(int ncol, int nlev, int nbnd, real3d const &bnd_flux_dn, real3d const &bnd_flux_up, real3d &bnd_flux_net) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(nbnd,nlev,ncol), YAKL_LAMBDA(int ibnd, int ilev, int icol) {
      bnd_flux_net(icol,ilev,ibnd) = bnd_flux_dn(icol,ilev,ibnd) - bnd_flux_up(icol,ilev,ibnd);
  });
}

#ifdef RRTMGP_ENABLE_KOKKOS
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2dk const &bnd_lims, real3dk const &spectral_flux, real3dk &byband_flux) {
  Kokkos::parallel_for( MDRangeP<3>({0,0,0}, {nbnd,nlev,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilev, int icol) {
    real bb_flux_s = 0.0;
    for (int igpt=bnd_lims(0,ibnd); igpt<=bnd_lims(1,ibnd); igpt++) {
      bb_flux_s += spectral_flux(icol,ilev,igpt);
    }
    byband_flux(icol,ilev,ibnd) = bb_flux_s;
  });
}
// Compute net flux
void net_byband(int ncol, int nlev, int nbnd, real3dk const &bnd_flux_dn, real3dk const &bnd_flux_up, real3dk &bnd_flux_net) {
  Kokkos::parallel_for( MDRangeP<3>({0,0,0}, {nbnd,nlev,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilev, int icol) {
      bnd_flux_net(icol,ilev,ibnd) = bnd_flux_dn(icol,ilev,ibnd) - bnd_flux_up(icol,ilev,ibnd);
  });
}
#endif
