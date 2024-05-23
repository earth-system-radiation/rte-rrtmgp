#pragma once
#include "rrtmgp_const.h"
#include "rrtmgp_conversion.h"

#ifdef RRTMGP_ENABLE_YAKL
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &bnd_lims, real3d const &spectral_flux, real3d &byband_flux);
void net_byband(int ncol, int nlev, int nbnd, real3d const &bnd_flux_dn, real3d const &bnd_flux_up, real3d &bnd_flux_net);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
// Spectral reduction over all points
template <typename BndLimsT, typename SpectralT, typename BybandT>
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, BndLimsT const &bnd_lims,
                SpectralT const &spectral_flux, BybandT &byband_flux) {
  using RealT = typename SpectralT::non_const_value_type;
  Kokkos::parallel_for( conv::get_mdrp<3>({nbnd,nlev,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilev, int icol) {
    RealT bb_flux_s = 0.0;
    for (int igpt=bnd_lims(0,ibnd); igpt<=bnd_lims(1,ibnd); igpt++) {
      bb_flux_s += spectral_flux(icol,ilev,igpt);
    }
    byband_flux(icol,ilev,ibnd) = bb_flux_s;
  });
}
// Compute net flux
template <typename FluxDnT, typename FluxUpT, typename FluxNetT>
void net_byband(int ncol, int nlev, int nbnd,
                FluxDnT const &bnd_flux_dn, FluxUpT const &bnd_flux_up, FluxNetT &bnd_flux_net) {
  Kokkos::parallel_for( conv::get_mdrp<3>({nbnd,nlev,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilev, int icol) {
      bnd_flux_net(icol,ilev,ibnd) = bnd_flux_dn(icol,ilev,ibnd) - bnd_flux_up(icol,ilev,ibnd);
  });
}
#endif
