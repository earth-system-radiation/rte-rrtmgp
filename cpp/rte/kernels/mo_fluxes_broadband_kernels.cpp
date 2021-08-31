#include "mo_fluxes_broadband_kernels.h"

// Spectral reduction over all points
void sum_broadband(int ncol, int nlev, int ngpt, real3d const &spectral_flux, real2d &broadband_flux) {
  // do ilev = 1, nlev
  //   do icol = 1, ncol
  parallel_for( Bounds<2>(nlev,ncol) , YAKL_LAMBDA (int ilev, int icol) {
    real bb_flux_s = 0.0_wp;
    for (int igpt=1; igpt<=ngpt; igpt++) {
      bb_flux_s += spectral_flux(icol, ilev, igpt);
    }
    broadband_flux(icol, ilev) = bb_flux_s;
  });
}

// Net flux: Spectral reduction over all points
void net_broadband(int ncol, int nlev, int ngpt, real3d const &spectral_flux_dn, real3d const &spectral_flux_up, real2d &broadband_flux_net) {
  // do ilev = 1, nlev
  //   do icol = 1, ncol
  parallel_for( Bounds<2>(nlev,ncol) , YAKL_LAMBDA (int ilev, int icol) {
    real diff = spectral_flux_dn(icol, ilev, 1) - spectral_flux_up(icol, ilev, 1);
    broadband_flux_net(icol, ilev) = diff;
  });

  // do igpt = 2, ngpt
  //   do ilev = 1, nlev
  //     do icol = 1, ncol
  parallel_for( Bounds<3>({2,ngpt},nlev,ncol) , YAKL_DEVICE_LAMBDA (int igpt, int ilev, int icol) {
    real diff = spectral_flux_dn(icol, ilev, igpt) - spectral_flux_up(icol, ilev, igpt);
    yakl::atomicAdd( broadband_flux_net(icol,ilev) , diff );
  });
#ifdef RRTMGP_DEBUG
  std::cout << "WARNING: THIS ISN'T TESTED!\n";
  std::cout << __FILE__ << ": " << __LINE__ << std::endl;
#endif
}

// Net flux when bradband flux up and down are already available
void net_broadband(int ncol, int nlev, real2d const &flux_dn, real2d const &flux_up, real2d &broadband_flux_net) {
  // do ilev = 1, nlev
  //   do icol = 1, ncol
  parallel_for( Bounds<2>(nlev,ncol) , YAKL_LAMBDA (int ilev, int icol) {
     broadband_flux_net(icol,ilev) = flux_dn(icol,ilev) - flux_up(icol,ilev);
  });
#ifdef RRTMGP_DEBUG
  std::cout << "WARNING: THIS ISN'T TESTED!\n";
  std::cout << __FILE__ << ": " << __LINE__ << std::endl;
#endif
}
