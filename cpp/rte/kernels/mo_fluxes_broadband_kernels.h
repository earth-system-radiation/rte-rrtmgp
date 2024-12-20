#pragma once

#include "rrtmgp_const.h"
#include "rrtmgp_conversion.h"

#ifdef RRTMGP_ENABLE_YAKL
// Spectral reduction over all points
void sum_broadband(int ncol, int nlev, int ngpt, real3d const &spectral_flux, real2d const &broadband_flux);

// Net flux: Spectral reduction over all points
void net_broadband(int ncol, int nlev, int ngpt, real3d const &spectral_flux_dn, real3d const &spectral_flux_up, real2d const &broadband_flux_net);

// Net flux when bradband flux up and down are already available
void net_broadband(int ncol, int nlev, real2d const &flux_dn, real2d const &flux_up, real2d const &broadband_flux_net);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
// Spectral reduction over all points
template <typename SpectralT, typename BroadbandT>
void sum_broadband(int ncol, int nlev, int ngpt, SpectralT const &spectral_flux, BroadbandT const &broadband_flux) {
  using RealT = typename SpectralT::non_const_value_type;
  using DeviceT = typename SpectralT::device_type;
  using LayoutT = typename SpectralT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  // do ilev = 1, nlev
  //   do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ncol, nlev}) , KOKKOS_LAMBDA (int icol, int ilev) {
    RealT bb_flux_s = 0.0;
    for (int igpt=0; igpt<ngpt; igpt++) {
      bb_flux_s += spectral_flux(icol, ilev, igpt);
    }
    broadband_flux(icol, ilev) = bb_flux_s;
  }));
}

// Net flux: Spectral reduction over all points
template <typename SpectralDnT, typename SpectralUpT, typename BroadbandT>
void net_broadband(int ncol, int nlev, int ngpt, SpectralDnT const &spectral_flux_dn, SpectralUpT const &spectral_flux_up, BroadbandT const &broadband_flux_net) {
  using RealT = typename SpectralDnT::non_const_value_type;
  using DeviceT = typename SpectralDnT::device_type;
  using LayoutT = typename SpectralDnT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  // do ilev = 1, nlev
  //   do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ncol, nlev}) , KOKKOS_LAMBDA (int icol, int ilev) {
    RealT diff = spectral_flux_dn(icol, ilev, 0) - spectral_flux_up(icol, ilev, 0);
    broadband_flux_net(icol, ilev) = diff;
  }));

  // do igpt = 2, ngpt
  //   do ilev = 1, nlev
  //     do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ncol, nlev}) , KOKKOS_LAMBDA (int icol, int ilev) {
    for (int igpt=1; igpt<ngpt; igpt++) {
      RealT diff = spectral_flux_dn(icol, ilev, igpt) - spectral_flux_up(icol, ilev, igpt);
      broadband_flux_net(icol,ilev) += diff;
    }
  }));
#ifdef RRTMGP_DEBUG
  std::cout << "WARNING: THIS ISN'T TESTED!\n";
  std::cout << __FILE__ << ": " << __LINE__ << std::endl;
#endif
}

// Net flux when bradband flux up and down are already available
template <typename FluxDnT, typename FluxUpT, typename BroadbandT>
void net_broadband(int ncol, int nlev, FluxDnT const &flux_dn, FluxUpT const &flux_up, BroadbandT const &broadband_flux_net) {
  using DeviceT = typename FluxDnT::device_type;
  using LayoutT = typename FluxDnT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  // do ilev = 1, nlev
  //   do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ncol, nlev}) , KOKKOS_LAMBDA (int icol, int ilev) {
     broadband_flux_net(icol,ilev) = flux_dn(icol,ilev) - flux_up(icol,ilev);
  }));
}
#endif
