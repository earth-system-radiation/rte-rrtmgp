#include "mo_fluxes_byband_kernels.h"

// Spectral reduction over all points
void sum_byband(int ncol, int nlev, int ngpt, int nbnd, int2d const &bnd_lims, real3d const &spectral_flux, real3d &byband_flux) {
  parallel_for( Bounds<3>(nbnd,nlev,ncol) , YAKL_LAMBDA (int ibnd, int ilev, int icol) {
    real bb_flux_s = 0.0_wp;
    for (int igpt=bnd_lims(1,ibnd); igpt<=bnd_lims(2,ibnd); igpt++) {
      bb_flux_s += spectral_flux(icol,ilev,igpt);
    }
    byband_flux(icol,ilev,ibnd) = bb_flux_s;
  });
}
