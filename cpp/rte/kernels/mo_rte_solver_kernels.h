
#pragma once

#include "rrtmgp_const.h"



//   Lower-level longwave kernels
// ---------------------------------------------------------------
// Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
// See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
// This routine implements point-wise stencil, and has to be called in a loop
// ---------------------------------------------------------------
YAKL_INLINE void lw_source_noscat_stencil(int ncol, int nlay, int ngpt, int icol, int ilay, int igpt,
                                          real3d const &lay_source, real3d const &lev_source_up, real3d const &lev_source_dn,
                                          real3d const &tau, real3d const &trans, real3d const &source_dn, real3d const &source_up, real tau_thresh) {
  // Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
  //   is of order epsilon (smallest difference from 1. in working precision)
  //   Thanks to Peter Blossey
  const auto term1 = trans(icol,ilay,igpt);
  const auto term2 = tau(icol,ilay,igpt);
  real fact;
  if (term2 > tau_thresh) {
    fact = (1._wp - term1) / term2 - term1;
  } else {
    fact = term2 * ( 0.5_wp - 1._wp/3._wp*term2 );
  }
  // Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
  source_dn(icol,ilay,igpt) = (1._wp - term1) * lev_source_dn(icol,ilay,igpt) +
                              2._wp * fact * (lay_source(icol,ilay,igpt) - lev_source_dn(icol,ilay,igpt));
  source_up(icol,ilay,igpt) = (1._wp - term1) * lev_source_up(icol,ilay,igpt) +
                              2._wp * fact * (lay_source(icol,ilay,igpt) - lev_source_up(icol,ilay,igpt));
}



// Direct beam source for diffuse radiation in layers and at surface;
//   report direct beam as a byproduct
inline void sw_source_2str(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &Rdir, real3d const &Tdir,
                           real3d const &Tnoscat, real2d const &sfc_albedo, real3d const &source_up, real3d const &source_dn,
                           real2d const &source_sfc, real3d const &flux_dn_dir) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( KERNEL_NAME() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      for (int ilev=1; ilev<=nlay; ilev++) {
        source_up(icol,ilev,igpt)     =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        source_dn(icol,ilev,igpt)     =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        flux_dn_dir(icol,ilev+1,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        if (ilev == nlay) {
          source_sfc(icol,igpt) = flux_dn_dir(icol,nlay+1,igpt)*sfc_albedo(icol,igpt);
        }
      }
    });
  } else {
    // layer index = level index
    // previous level is up (+1)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( KERNEL_NAME() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      for (int ilev=nlay; ilev>=1; ilev--) {
        source_up(icol,ilev,igpt)   =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        source_dn(icol,ilev,igpt)   =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        flux_dn_dir(icol,ilev,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        if (ilev ==    1) {
          source_sfc(icol,igpt) = flux_dn_dir(icol,    1,igpt)*sfc_albedo(icol,igpt);
        }
      }
    });
  }
}



// Longwave no-scattering transport
inline void lw_transport_noscat(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real3d const &trans,
                                real2d const &sfc_albedo, real3d const &source_dn, real3d const &source_up, real2d const &source_sfc, 
                                real3d const &radn_up, real3d const &radn_dn) {
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  if (top_at_1) {
    // Top of domain is index 1
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( KERNEL_NAME() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      // Downward propagation
      for (int ilev=2; ilev<=nlay+1; ilev++) {
        radn_dn(icol,ilev,igpt) = trans(icol,ilev-1,igpt)*radn_dn(icol,ilev-1,igpt) + source_dn(icol,ilev-1,igpt);
      }

      // Surface reflection and emission
      radn_up(icol,nlay+1,igpt) = radn_dn(icol,nlay+1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt);

      // Upward propagation
      for (int ilev=nlay; ilev>=1; ilev--) {
        radn_up(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_up(icol,ilev+1,igpt) + source_up(icol,ilev,igpt);
      }
    });
  } else {
    // Top of domain is index nlay+1
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( KERNEL_NAME() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      // Downward propagation
      for (int ilev=nlay; ilev>=1; ilev--) {
        radn_dn(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_dn(icol,ilev+1,igpt) + source_dn(icol,ilev,igpt);
      }

      // Surface reflection and emission
      radn_up(icol,     1,igpt) = radn_dn(icol,     1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt);

      // Upward propagation
      for (int ilev=2; ilev<=nlay+1; ilev++) {
        radn_up(icol,ilev,igpt) = trans(icol,ilev-1,igpt) * radn_up(icol,ilev-1,igpt) +  source_up(icol,ilev-1,igpt);
      }
    });
  }
}



//   Lower-level shortwave kernels
//
// Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
//    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
//
// Equations are developed in Meador and Weaver, 1980,
//    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
inline void sw_two_stream(int ncol, int nlay, int ngpt, real1d const &mu0, real3d const &tau,
                          real3d const &w0, real3d const &g, real3d &Rdif, real3d const &Tdif,
                          real3d const &Rdir, real3d const &Tdir, real3d const &Tnoscat) {
  using yakl::intrinsics::merge;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real1d mu0_inv("mu0_inv",ncol);

  real eps = std::numeric_limits<real>::epsilon();

  parallel_for( KERNEL_NAME() , ncol , YAKL_LAMBDA (int icol) {
    mu0_inv(icol) = 1._wp/mu0(icol);
  });

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  parallel_for( KERNEL_NAME() , SimpleBounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    // Zdunkowski Practical Improved Flux Method "PIFM"
    //  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    real gamma1= (8._wp - w0(icol,ilay,igpt) * (5._wp + 3._wp * g(icol,ilay,igpt))) * .25_wp;
    real gamma2=  3._wp *(w0(icol,ilay,igpt) * (1._wp -         g(icol,ilay,igpt))) * .25_wp;
    real gamma3= (2._wp - 3._wp * mu0(icol)  *                  g(icol,ilay,igpt) ) * .25_wp;
    real gamma4=  1._wp - gamma3;

    real alpha1 = gamma1 * gamma4 + gamma2 * gamma3;           // Eq. 16
    real alpha2 = gamma1 * gamma3 + gamma2 * gamma4;           // Eq. 17
    // Written to encourage vectorization of exponential, square root
    // Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
    //   k = 0 for isotropic, conservative scattering; this lower limit on k
    //   gives relative error with respect to conservative solution
    //   of < 0.1% in Rdif down to tau = 10^-9
    real k = sqrt(max((gamma1 - gamma2) * 
                      (gamma1 + gamma2),  
                      1.e-12_wp));
    real exp_minusktau = exp(-tau(icol,ilay,igpt)*k);

    // Diffuse reflection and transmission
    real exp_minus2ktau = exp_minusktau * exp_minusktau;

    // Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    real RT_term = 1._wp / (k      * (1._wp + exp_minus2ktau)  + 
                            gamma1 * (1._wp - exp_minus2ktau) );

    // Equation 25
    Rdif(icol,ilay,igpt) = RT_term * gamma2 * (1._wp - exp_minus2ktau);

    // Equation 26
    Tdif(icol,ilay,igpt) = RT_term * 2._wp * k * exp_minusktau;

    // Transmittance of direct, unscattered beam. Also used below
    Tnoscat(icol,ilay,igpt) = exp(-tau(icol,ilay,igpt)*mu0_inv(icol));

    // Direct reflect and transmission
    real k_mu     = k * mu0(icol);
    real k_gamma3 = k * gamma3;
    real k_gamma4 = k * gamma4;

    // Equation 14, multiplying top and bottom by exp(-k*tau)
    //   and rearranging to avoid div by 0.
    RT_term =  w0(icol,ilay,igpt) * RT_term/merge(1._wp - k_mu*k_mu, 
                                                  eps,    
                                                  abs(1._wp - k_mu*k_mu) >= eps);

    Rdir(icol,ilay,igpt) = RT_term  *                                    
       ((1._wp - k_mu) * (alpha2 + k_gamma3)                  - 
        (1._wp + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau - 
        2.0_wp * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau  * Tnoscat(icol,ilay,igpt));

    // Equation 15, multiplying top and bottom by exp(-k*tau),
    //   multiplying through by exp(-tau/mu0) to
    //   prefer underflow to overflow
    // Omitting direct transmittance
    Tdir(icol,ilay,igpt) = 
             -RT_term * ((1._wp + k_mu) * (alpha1 + k_gamma4) * Tnoscat(icol,ilay,igpt) - 
                         (1._wp - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * Tnoscat(icol,ilay,igpt) - 
                          2.0_wp * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau );

  });
}



void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &flux_dn);



void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &inc_flux, real1d const &factor, real3d const &flux_dn);


// Upper boundary condition
void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &inc_flux, real3d const &flux_dn);



// Transport of diffuse radiation through a vertically layered atmosphere.
//   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
//   This routine is shared by longwave and shortwave
void adding(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &albedo_sfc, real3d const &rdif, real3d const &tdif,
            real3d const &src_dn, real3d const &src_up, real2d const &src_sfc, real3d const &flux_up, real3d const &flux_dn);




void sw_solver_2stream(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real3d const &ssa, real3d const &g,
                       real1d const &mu0, real2d const &sfc_alb_dir, real2d const &sfc_alb_dif, real3d const &flux_up,
                       real3d const &flux_dn, real3d const &flux_dir);



// Top-level longwave kernels
//
// LW fluxes, no scattering, mu (cosine of integration angle) specified by column
//   Does radiation calculation at user-supplied angles; converts radiances to flux
//   using user-supplied weights
void lw_solver_noscat(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &D, real1d const &weights, int weight_ind, real3d const &tau,
                      real3d const &lay_source, real3d const &lev_source_inc, real3d const &lev_source_dec,
                      real2d const &sfc_emis, real2d const &sfc_src, real3d const &radn_up, real3d const &radn_dn);



// LW transport, no scattering, multi-angle quadrature
//   Users provide a set of weights and quadrature angles
//   Routine sums over single-angle solutions for each sets of angles/weights
void lw_solver_noscat_GaussQuad(int ncol, int nlay, int ngpt, bool top_at_1, int nmus, real1d const &Ds, real1d const &weights, 
                                real3d const &tau, real3d const &lay_source, real3d const &lev_source_inc, real3d const &lev_source_dec,
                                real2d const &sfc_emis, real2d const &sfc_src, real3d const &flux_up, real3d const &flux_dn);



// Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
//   This version straight from ECRAD
//   Source is provided as W/m2-str; factor of pi converts to flux units
void lw_source_2str(int ncol, int nlay, int ngpt, bool top_at_1, real2d const &sfc_emis, real2d const &sfc_src,
                    real3d const &lay_source, real3d const &lev_source, real3d const &gamma1, real3d const &gamma2,
                    real3d const &rdif, real3d const &tdif, real3d const &tau, real3d const &source_dn, real3d const &source_up,
                    real2d const &source_sfc);



// Source function combination
// RRTMGP provides two source functions at each level
//   using the spectral mapping from each of the adjascent layers.
//   Need to combine these for use in two-stream calculation.
void lw_combine_sources(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &lev_src_inc, real3d const &lev_src_dec,
                        real3d const &lev_source);



// Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
//    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
// Equations are developed in Meador and Weaver, 1980,
//    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
void lw_two_stream(int ncol, int nlay, int ngpt, real3d const &tau, real3d const &w0, real3d const &g, real3d const &gamma1,
                   real3d const &gamma2, real3d const &Rdif, real3d const &Tdif);



//   Top-level shortwave kernels
// -------------------------------------------------------------------------------------------------
//   Extinction-only i.e. solar direct beam
void sw_solver_noscat(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real1d const &mu0, real3d const &flux_dir);



// Longwave two-stream calculation:
//   combine RRTMGP-specific sources at levels
//   compute layer reflectance, transmittance
//   compute total source function at levels using linear-in-tau
//   transport
void lw_solver_2stream(int ncol, int nlay, int ngpt, bool top_at_1, real3d const &tau, real3d const &ssa, real3d const &g,
                       real3d const &lay_source, real3d const &lev_source_inc, real3d const &lev_source_dec,
                       real2d const &sfc_emis, real2d const &sfc_src, real3d const &flux_up, real3d const &flux_dn);





