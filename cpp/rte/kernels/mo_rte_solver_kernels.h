
#pragma once

#include "rrtmgp_const.h"
#include "rrtmgp_conversion.h"

#ifdef RRTMGP_ENABLE_YAKL
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
    fact = (1. - term1) / term2 - term1;
  } else {
    fact = term2 * ( 0.5 - 1./3.*term2 );
  }
  // Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
  source_dn(icol,ilay,igpt) = (1. - term1) * lev_source_dn(icol,ilay,igpt) +
                              2. * fact * (lay_source(icol,ilay,igpt) - lev_source_dn(icol,ilay,igpt));
  source_up(icol,ilay,igpt) = (1. - term1) * lev_source_up(icol,ilay,igpt) +
                              2. * fact * (lay_source(icol,ilay,igpt) - lev_source_up(icol,ilay,igpt));
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
    TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      for (int ilev=1; ilev<=nlay; ilev++) {
        source_up(icol,ilev,igpt)     =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        source_dn(icol,ilev,igpt)     =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        flux_dn_dir(icol,ilev+1,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        if (ilev == nlay) {
          source_sfc(icol,igpt) = flux_dn_dir(icol,nlay+1,igpt)*sfc_albedo(icol,igpt);
        }
      }
    }));
  } else {
    #ifdef RRTMGP_CPU_KERNELS
      #ifdef YAKL_AUTO_PROFILE
        auto timername = std::string(YAKL_AUTO_LABEL());
        yakl::timer_start(timername.c_str());
      #endif
      #ifdef YAKL_ARCH_OPENMP
        #pragma omp parallel for
      #endif
      for (int igpt = 1; igpt <= ngpt; igpt++) {
        for (int ilev=nlay; ilev>=1; ilev--) {
          for (int icol = 1; icol <= ncol; icol++) {
            source_up(icol,ilev,igpt)   =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
            source_dn(icol,ilev,igpt)   =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
            flux_dn_dir(icol,ilev,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
            if (ilev ==    1) {
              source_sfc(icol,igpt) = flux_dn_dir(icol,    1,igpt)*sfc_albedo(icol,igpt);
            }
          }
        }
      }
      #ifdef YAKL_AUTO_PROFILE
        yakl::timer_stop(timername.c_str());
      #endif
    #else
      // layer index = level index
      // previous level is up (+1)
      // do igpt = 1, ngpt
      //   do icol = 1, ncol
      TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
        for (int ilev=nlay; ilev>=1; ilev--) {
          source_up(icol,ilev,igpt)   =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
          source_dn(icol,ilev,igpt)   =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
          flux_dn_dir(icol,ilev,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
          if (ilev ==    1) {
            source_sfc(icol,igpt) = flux_dn_dir(icol,    1,igpt)*sfc_albedo(icol,igpt);
          }
        }
      }));
    #endif
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
    TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
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
    }));
  } else {
    // Top of domain is index nlay+1
    #ifdef RRTMGP_CPU_KERNELS
      #ifdef YAKL_AUTO_PROFILE
        auto timername = std::string(YAKL_AUTO_LABEL());
        yakl::timer_start(timername.c_str());
      #endif
      #ifdef YAKL_ARCH_OPENMP
        #pragma omp parallel for
      #endif
      for (int igpt = 1; igpt <= ngpt; igpt++) {
        // Downward propagation
        for (int ilev=nlay; ilev>=1; ilev--) {
          for (int icol = 1; icol <= ncol; icol++) {
            radn_dn(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_dn(icol,ilev+1,igpt) + source_dn(icol,ilev,igpt);
          }
        }
        // Surface reflection and emission
        for (int icol = 1; icol <= ncol; icol++) {
          radn_up(icol,     1,igpt) = radn_dn(icol,     1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt);
        }
        // Upward propagation
        for (int ilev=2; ilev<=nlay+1; ilev++) {
          for (int icol = 1; icol <= ncol; icol++) {
            radn_up(icol,ilev,igpt) = trans(icol,ilev-1,igpt) * radn_up(icol,ilev-1,igpt) +  source_up(icol,ilev-1,igpt);
          }
        }
      }
      #ifdef YAKL_AUTO_PROFILE
        yakl::timer_stop(timername.c_str());
      #endif
    #else
      // do igpt = 1, ngpt
      //   do icol = 1, ncol
      TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
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
      }));
    #endif
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

  TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , ncol , YAKL_LAMBDA (int icol) {
    mu0_inv(icol) = 1./mu0(icol);
  }));

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  TIMED_KERNEL(parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngpt,nlay,ncol) , YAKL_LAMBDA (int igpt, int ilay, int icol) {
    // Zdunkowski Practical Improved Flux Method "PIFM"
    //  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    real gamma1= (8. - w0(icol,ilay,igpt) * (5. + 3. * g(icol,ilay,igpt))) * .25;
    real gamma2=  3. *(w0(icol,ilay,igpt) * (1. -         g(icol,ilay,igpt))) * .25;
    real gamma3= (2. - 3. * mu0(icol)  *                  g(icol,ilay,igpt) ) * .25;
    real gamma4=  1. - gamma3;

    real alpha1 = gamma1 * gamma4 + gamma2 * gamma3;           // Eq. 16
    real alpha2 = gamma1 * gamma3 + gamma2 * gamma4;           // Eq. 17
    // Written to encourage vectorization of exponential, square root
    // Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
    //   k = 0 for isotropic, conservative scattering; this lower limit on k
    //   gives relative error with respect to conservative solution
    //   of < 0.1% in Rdif down to tau = 10^-9
    real k = sqrt(std::max((gamma1 - gamma2) *
                           (gamma1 + gamma2),
                           1.e-12));
    real exp_minusktau = exp(-tau(icol,ilay,igpt)*k);

    // Diffuse reflection and transmission
    real exp_minus2ktau = exp_minusktau * exp_minusktau;

    // Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    real RT_term = 1. / (k      * (1. + exp_minus2ktau)  +
                            gamma1 * (1. - exp_minus2ktau) );

    // Equation 25
    Rdif(icol,ilay,igpt) = RT_term * gamma2 * (1. - exp_minus2ktau);

    // Equation 26
    Tdif(icol,ilay,igpt) = RT_term * 2. * k * exp_minusktau;

    // Transmittance of direct, unscattered beam. Also used below
    Tnoscat(icol,ilay,igpt) = exp(-tau(icol,ilay,igpt)*mu0_inv(icol));

    // Direct reflect and transmission
    real k_mu     = k * mu0(icol);
    real k_gamma3 = k * gamma3;
    real k_gamma4 = k * gamma4;

    // Equation 14, multiplying top and bottom by exp(-k*tau)
    //   and rearranging to avoid div by 0.
    RT_term =  w0(icol,ilay,igpt) * RT_term/merge(1. - k_mu*k_mu,
                                                  eps,
                                                  std::abs(1. - k_mu*k_mu) >= eps);

    Rdir(icol,ilay,igpt) = RT_term  *
       ((1. - k_mu) * (alpha2 + k_gamma3)                  -
        (1. + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau -
        2.0 * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau  * Tnoscat(icol,ilay,igpt));

    // Equation 15, multiplying top and bottom by exp(-k*tau),
    //   multiplying through by exp(-tau/mu0) to
    //   prefer underflow to overflow
    // Omitting direct transmittance
    Tdir(icol,ilay,igpt) =
             -RT_term * ((1. + k_mu) * (alpha1 + k_gamma4) * Tnoscat(icol,ilay,igpt) -
                         (1. - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * Tnoscat(icol,ilay,igpt) -
                          2.0 * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau );

  }));
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
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
//   Lower-level longwave kernels
// ---------------------------------------------------------------
// Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
// See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
// This routine implements point-wise stencil, and has to be called in a loop
// ---------------------------------------------------------------
template <typename LaySourceT, typename LevUpT, typename LevDnT, typename TauT, typename TransT,
          typename SourceDnT, typename SourceUpT>
KOKKOS_INLINE_FUNCTION
void lw_source_noscat_stencil(
  int ncol, int nlay, int ngpt, int icol, int ilay, int igpt,
  LaySourceT const &lay_source, LevUpT const &lev_source_up, LevDnT const &lev_source_dn,
  TauT const &tau, TransT const &trans, SourceDnT const &source_dn, SourceUpT const &source_up,
  typename TauT::non_const_value_type tau_thresh)
{
  using RealT = typename TauT::non_const_value_type;
  // Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
  //   is of order epsilon (smallest difference from 1. in working precision)
  //   Thanks to Peter Blossey
  const auto term1 = trans(icol,ilay,igpt);
  const auto term2 = tau(icol,ilay,igpt);
  RealT fact;
  if (term2 > tau_thresh) {
    fact = (1. - term1) / term2 - term1;
  } else {
    fact = term2 * ( 0.5 - 1./3.*term2 );
  }
  // Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
  source_dn(icol,ilay,igpt) = (1. - term1) * lev_source_dn(icol,ilay,igpt) +
                              2. * fact * (lay_source(icol,ilay,igpt) - lev_source_dn(icol,ilay,igpt));
  source_up(icol,ilay,igpt) = (1. - term1) * lev_source_up(icol,ilay,igpt) +
                              2. * fact * (lay_source(icol,ilay,igpt) - lev_source_up(icol,ilay,igpt));
}

// Direct beam source for diffuse radiation in layers and at surface;
//   report direct beam as a byproduct
template <typename RdirT, typename TdirT, typename TnoscatT, typename SfcAlbedoT, typename SourceUpT,
          typename SourceDnT, typename SourceSfcT, typename FluxDnDirT>
inline void sw_source_2str(int ncol, int nlay, int ngpt, bool top_at_1, RdirT const &Rdir, TdirT const &Tdir,
                           TnoscatT const &Tnoscat, SfcAlbedoT const &sfc_albedo, SourceUpT const &source_up,
                           SourceDnT const &source_dn, SourceSfcT const &source_sfc, FluxDnDirT const &flux_dn_dir) {
  using DeviceT = typename RdirT::device_type;
  using LayoutT = typename RdirT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      for (int ilev=0; ilev<nlay; ilev++) {
        source_up(icol,ilev,igpt)     =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        source_dn(icol,ilev,igpt)     =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        flux_dn_dir(icol,ilev+1,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt);
        if (ilev == nlay-1) {
          source_sfc(icol,igpt) = flux_dn_dir(icol,nlay,igpt)*sfc_albedo(icol,igpt);
        }
      }
    }));
  } else {
    // layer index = level index
    // previous level is up (+1)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      for (int ilev=nlay-1; ilev>=0; ilev--) {
        source_up(icol,ilev,igpt)   =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        source_dn(icol,ilev,igpt)   =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        flux_dn_dir(icol,ilev,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt);
        if (ilev ==    0) {
          source_sfc(icol,igpt) = flux_dn_dir(icol,    0,igpt)*sfc_albedo(icol,igpt);
        }
      }
    }));
  }
}

// Longwave no-scattering transport
template <typename TauT, typename TransT, typename SfcAlbedoT, typename SourceDnT, typename SourceUpT,
          typename SourceSfcT, typename RadnUpT, typename RadnDnT>
inline void lw_transport_noscat(int ncol, int nlay, int ngpt, bool top_at_1, TauT const &tau, TransT const &trans,
                                SfcAlbedoT const &sfc_albedo, SourceDnT const &source_dn, SourceUpT const &source_up,
                                SourceSfcT const &source_sfc, RadnUpT const &radn_up, RadnDnT const &radn_dn) {
  using DeviceT = typename TauT::device_type;
  using LayoutT = typename TauT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  if (top_at_1) {
    // Top of domain is index 1
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      // Downward propagation
      for (int ilev=1; ilev<nlay+1; ilev++) {
        radn_dn(icol,ilev,igpt) = trans(icol,ilev-1,igpt)*radn_dn(icol,ilev-1,igpt) + source_dn(icol,ilev-1,igpt);
      }

      // Surface reflection and emission
      radn_up(icol,nlay,igpt) = radn_dn(icol,nlay,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt);

      // Upward propagation
      for (int ilev=nlay-1; ilev>=0; ilev--) {
        radn_up(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_up(icol,ilev+1,igpt) + source_up(icol,ilev,igpt);
      }
    }));
  } else {
    // Top of domain is index nlay+1
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    Kokkos::Array<int, 2> dims2_ncol_ngpt = {ncol,ngpt};
    const int dims2_tot = ncol*ngpt;
    TIMED_KERNEL(Kokkos::parallel_for( dims2_tot , KOKKOS_LAMBDA (int idx) {
      int icol, igpt;
      conv::unflatten_idx_left(idx, dims2_ncol_ngpt, icol, igpt);
      // Downward propagation
      for (int ilev=nlay-1; ilev>=0; ilev--) {
        radn_dn(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_dn(icol,ilev+1,igpt) + source_dn(icol,ilev,igpt);
      }

      // Surface reflection and emission
      radn_up(icol,0,igpt) = radn_dn(icol,0,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt);

      // Upward propagation
      for (int ilev=1; ilev<nlay+1; ilev++) {
        radn_up(icol,ilev,igpt) = trans(icol,ilev-1,igpt) * radn_up(icol,ilev-1,igpt) +  source_up(icol,ilev-1,igpt);
      }
    }));
  }
}

//   Lower-level shortwave kernels
//
// Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
//    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
//
// Equations are developed in Meador and Weaver, 1980,
//    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
template <typename Mu0T, typename TauT, typename W0T, typename GT, typename RdifT, typename TdifT,
          typename RdirT, typename TdirT, typename TnoscatT>
inline void sw_two_stream(int ncol, int nlay, int ngpt, Mu0T const &mu0, TauT const &tau,
                          W0T const &w0, GT const &g, RdifT &Rdif, TdifT const &Tdif,
                          RdirT const &Rdir, TdirT const &Tdir, TnoscatT const &Tnoscat) {
  using conv::merge;
  using RealT = typename TauT::non_const_value_type;
  using LayoutT = typename TauT::array_layout;
  using DeviceT = typename TauT::device_type;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;
  using pool = conv::MemPoolSingleton<RealT, DeviceT>;

  auto mu0_inv = pool::template alloc<RealT>(ncol);

  RealT eps = std::numeric_limits<RealT>::epsilon();

  TIMED_KERNEL(Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
    mu0_inv(icol) = 1./mu0(icol);
  }));

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ngpt,nlay,ncol}) , KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    // Zdunkowski Practical Improved Flux Method "PIFM"
    //  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
    RealT gamma1= (8. - w0(icol,ilay,igpt) * (5. + 3. * g(icol,ilay,igpt))) * .25;
    RealT gamma2=  3. *(w0(icol,ilay,igpt) * (1. -      g(icol,ilay,igpt))) * .25;
    RealT gamma3= (2. - 3. * mu0(icol)  *               g(icol,ilay,igpt) ) * .25;
    RealT gamma4=  1. - gamma3;

    RealT alpha1 = gamma1 * gamma4 + gamma2 * gamma3;           // Eq. 16
    RealT alpha2 = gamma1 * gamma3 + gamma2 * gamma4;           // Eq. 17
    // Written to encourage vectorization of exponential, square root
    // Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
    //   k = 0 for isotropic, conservative scattering; this lower limit on k
    //   gives relative error with respect to conservative solution
    //   of < 0.1% in Rdif down to tau = 10^-9
    RealT k = sqrt(Kokkos::fmax((gamma1 - gamma2) *
                           (gamma1 + gamma2),
                           1.e-12));
    RealT exp_minusktau = exp(-tau(icol,ilay,igpt)*k);

    // Diffuse reflection and transmission
    RealT exp_minus2ktau = exp_minusktau * exp_minusktau;

    // Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
    RealT RT_term = 1. / (k      * (1. + exp_minus2ktau)  +
                            gamma1 * (1. - exp_minus2ktau) );

    // Equation 25
    Rdif(icol,ilay,igpt) = RT_term * gamma2 * (1. - exp_minus2ktau);

    // Equation 26
    Tdif(icol,ilay,igpt) = RT_term * 2. * k * exp_minusktau;

    // Transmittance of direct, unscattered beam. Also used below
    Tnoscat(icol,ilay,igpt) = exp(-tau(icol,ilay,igpt)*mu0_inv(icol));

    // Direct reflect and transmission
    RealT k_mu     = k * mu0(icol);
    RealT k_gamma3 = k * gamma3;
    RealT k_gamma4 = k * gamma4;

    // Equation 14, multiplying top and bottom by exp(-k*tau)
    //   and rearranging to avoid div by 0.
    RT_term =  w0(icol,ilay,igpt) * RT_term/merge(1. - k_mu*k_mu,
                                                  eps,
                                                  Kokkos::fabs(1. - k_mu*k_mu) >= eps);

    Rdir(icol,ilay,igpt) = RT_term  *
       ((1. - k_mu) * (alpha2 + k_gamma3)                  -
        (1. + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau -
        2.0 * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau  * Tnoscat(icol,ilay,igpt));

    // Equation 15, multiplying top and bottom by exp(-k*tau),
    //   multiplying through by exp(-tau/mu0) to
    //   prefer underflow to overflow
    // Omitting direct transmittance
    Tdir(icol,ilay,igpt) =
             -RT_term * ((1. + k_mu) * (alpha1 + k_gamma4) * Tnoscat(icol,ilay,igpt) -
                         (1. - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * Tnoscat(icol,ilay,igpt) -
                          2.0 * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau );

  }));

  pool::dealloc(mu0_inv.data(), mu0_inv.size());
}

template <typename FluxDnT>
void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, FluxDnT const &flux_dn) {
  using DeviceT = typename FluxDnT::device_type;
  using LayoutT = typename FluxDnT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  //   Upper boundary condition
  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      flux_dn(icol,      0, igpt)  = 0;
    }));
  } else {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      flux_dn(icol, nlay, igpt)  = 0;
    }));
  }
}

template <typename IncFluxT, typename FactorT, typename FluxDnT>
void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, IncFluxT const &inc_flux, FactorT const &factor, FluxDnT const &flux_dn) {
  using DeviceT = typename FluxDnT::device_type;
  using LayoutT = typename FluxDnT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      flux_dn(icol,      0, igpt)  = inc_flux(icol,igpt) * factor(icol);
    }));
  } else {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      flux_dn(icol, nlay, igpt)  = inc_flux(icol,igpt) * factor(icol);
    }));
  }
}


// Upper boundary condition
template <typename IncFluxT, typename FluxDnT>
void apply_BC(int ncol, int nlay, int ngpt, bool top_at_1, IncFluxT const &inc_flux, FluxDnT const &flux_dn) {
  using DeviceT = typename FluxDnT::device_type;
  using LayoutT = typename FluxDnT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  //   Upper boundary condition
  if (top_at_1) {
    //$acc  parallel loop collapse(2)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      flux_dn(icol,      0, igpt)  = inc_flux(icol,igpt);
    }));
  } else {
    //$acc  parallel loop collapse(2)
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      flux_dn(icol, nlay, igpt)  = inc_flux(icol,igpt);
    }));
  }
}

// Transport of diffuse radiation through a vertically layered atmosphere.
//   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
//   This routine is shared by longwave and shortwave
template <typename AlbedoSfcT, typename RdifT, typename TdifT, typename SrcDnT, typename SrcUpT,
          typename SrcSfcT, typename FluxUpT, typename FluxDnT>
void adding(int ncol, int nlay, int ngpt, bool top_at_1, AlbedoSfcT const &albedo_sfc, RdifT const &rdif, TdifT const &tdif,
            SrcDnT const &src_dn, SrcUpT const &src_up, SrcSfcT const &src_sfc, FluxUpT const &flux_up, FluxDnT const &flux_dn) {
  using RealT   = typename TdifT::non_const_value_type;
  using LayoutT = typename TdifT::array_layout;
  using DeviceT = typename TdifT::device_type;
  using ureal3d_t = conv::Unmanaged<Kokkos::View<RealT***, LayoutT, DeviceT>>;
  using pool    = conv::MemPoolSingleton<RealT, DeviceT>;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  const int dsize1 = ncol*(nlay+1)*ngpt;
  const int dsize2 = ncol*nlay*ngpt;
  RealT* data = pool::template alloc_raw<RealT>(dsize1*2 + dsize2), *dcurr = data;
  ureal3d_t albedo(dcurr,ncol,nlay+1,ngpt); dcurr += dsize1;
  ureal3d_t src   (dcurr,ncol,nlay+1,ngpt); dcurr += dsize1;
  ureal3d_t denom (dcurr,ncol,nlay  ,ngpt); dcurr += dsize2;

  // Indexing into arrays for upward and downward propagation depends on the vertical
  //   orientation of the arrays (whether the domain top is at the first or last index)
  // We write the loops out explicitly so compilers will have no trouble optimizing them.
  if (top_at_1) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      int ilev = nlay;
      // Albedo of lowest level is the surface albedo...
      albedo(icol,ilev,igpt)  = albedo_sfc(icol,igpt);
      // ... and source of diffuse radiation is surface emission
      src(icol,ilev,igpt) = src_sfc(icol,igpt);

      // From bottom to top of atmosphere --
      //   compute albedo and source of upward radiation
      for (ilev=nlay-1; ilev>=0; ilev--) {
        denom(icol,ilev,igpt) = 1./(1. - rdif(icol,ilev,igpt)*albedo(icol,ilev+1,igpt));    // Eq 10
        albedo(icol,ilev,igpt) = rdif(icol,ilev,igpt) +
                                 tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(icol,ilev+1,igpt) * denom(icol,ilev,igpt); // Equation 9
        // Equation 11 -- source is emitted upward radiation at top of layer plus
        //   radiation emitted at bottom of layer,
        //   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
        src(icol,ilev,igpt) =  src_up(icol, ilev, igpt) +
                               tdif(icol,ilev,igpt) * denom(icol,ilev,igpt) *
                               (src(icol,ilev+1,igpt) + albedo(icol,ilev+1,igpt)*src_dn(icol,ilev,igpt));
      }

      // Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = 0;
      flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // ... reflection of incident diffuse and
                                src(icol,ilev,igpt);                                  // emission from below

      // From the top of the atmosphere downward -- compute fluxes
      for (ilev = 1; ilev <= nlay; ilev++) {
        flux_dn(icol,ilev,igpt) = (tdif(icol,ilev-1,igpt)*flux_dn(icol,ilev-1,igpt) +   // Equation 13
                                  rdif(icol,ilev-1,igpt)*src(icol,ilev,igpt) +
                                  src_dn(icol,ilev-1,igpt)) * denom(icol,ilev-1,igpt);
        flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // Equation 12
                                  src(icol,ilev,igpt);
      }
    }));

  } else {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      int ilev = 0;
      // Albedo of lowest level is the surface albedo...
      albedo(icol,ilev,igpt)  = albedo_sfc(icol,igpt);
      // ... and source of diffuse radiation is surface emission
      src(icol,ilev,igpt) = src_sfc(icol,igpt);

      // From bottom to top of atmosphere --
      //   compute albedo and source of upward radiation
      for (ilev = 0; ilev < nlay; ilev++) {
        denom (icol,ilev  ,igpt) = 1./(1. - rdif(icol,ilev,igpt)*albedo(icol,ilev,igpt));                // Eq 10
        albedo(icol,ilev+1,igpt) = rdif(icol,ilev,igpt) +
          tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(icol,ilev,igpt) * denom(icol,ilev,igpt); // Equation 9
        // Equation 11 -- source is emitted upward radiation at top of layer plus
        //   radiation emitted at bottom of layer,
        //   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
        src(icol,ilev+1,igpt) =  src_up(icol, ilev, igpt) +
          tdif(icol,ilev,igpt) * denom(icol,ilev,igpt) *
          (src(icol,ilev,igpt) + albedo(icol,ilev,igpt)*src_dn(icol,ilev,igpt));
      }

      // Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = nlay;
      flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // ... reflection of incident diffuse and
        src(icol,ilev,igpt);                          // scattering by the direct beam below

      // From the top of the atmosphere downward -- compute fluxes
      for (ilev=nlay-1; ilev >= 0; ilev--) {
        flux_dn(icol,ilev,igpt) = (tdif(icol,ilev,igpt)*flux_dn(icol,ilev+1,igpt) +   // Equation 13
                                   rdif(icol,ilev,igpt)*src(icol,ilev,igpt) +
                                   src_dn(icol, ilev, igpt)) * denom(icol,ilev,igpt);
        flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(icol,ilev,igpt) +  // Equation 12
                                  src(icol,ilev,igpt);

      }
    }));
  }

  pool::dealloc(data, dcurr - data);
}

template <typename TauT, typename SsaT, typename GT, typename Mu0T, typename SfcDirT, typename SfcDifT,
          typename FluxUpT, typename FluxDnT, typename FluxDirT>
void sw_solver_2stream(int ncol, int nlay, int ngpt, bool top_at_1, TauT const &tau, SsaT const &ssa, GT const &g,
                       Mu0T const &mu0, SfcDirT const &sfc_alb_dir, SfcDifT const &sfc_alb_dif, FluxUpT const &flux_up,
                       FluxDnT const &flux_dn, FluxDirT const &flux_dir) {
  using RealT   = typename TauT::non_const_value_type;
  using LayoutT = typename TauT::array_layout;
  using DeviceT = typename TauT::device_type;
  using pool = conv::MemPoolSingleton<RealT, DeviceT>;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  using ureal2d_t = conv::Unmanaged<Kokkos::View<RealT**, LayoutT, DeviceT>>;
  using ureal3d_t = conv::Unmanaged<Kokkos::View<RealT***, LayoutT, DeviceT>>;

  const int dsize1 = ncol*nlay*ngpt;
  const int dsize2 = ncol*ngpt;
  RealT* data = pool::template alloc_raw<RealT>(dsize1*7 + dsize2), *dcurr = data;
  ureal3d_t Rdif      (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t Tdif      (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t Rdir      (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t Tdir      (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t Tnoscat   (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t source_up (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t source_dn (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal2d_t source_srf(dcurr,ncol     ,ngpt); dcurr += dsize2;

  // Cell properties: transmittance and reflectance for direct and diffuse radiation
  sw_two_stream(ncol, nlay, ngpt, mu0,
                tau , ssa , g   ,
                Rdif, Tdif, Rdir, Tdir, Tnoscat);

  sw_source_2str(ncol, nlay, ngpt, top_at_1,
                 Rdir, Tdir, Tnoscat, sfc_alb_dir,
                 source_up, source_dn, source_srf, flux_dir);

  adding(ncol, nlay, ngpt, top_at_1,
         sfc_alb_dif, Rdif, Tdif,
         source_dn, source_up, source_srf, flux_up, flux_dn);

  // adding computes only diffuse flux; flux_dn is total
  //
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay+1
  //     do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ngpt,nlay+1,ncol}) , KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    flux_dn(icol,ilay,igpt) = flux_dn(icol,ilay,igpt) + flux_dir(icol,ilay,igpt);
  }));

  pool::dealloc(data, dcurr - data);
}

// Top-level longwave kernels
//
// LW fluxes, no scattering, mu (cosine of integration angle) specified by column
//   Does radiation calculation at user-supplied angles; converts radiances to flux
//   using user-supplied weights
template <typename DT, typename WeightsT, typename TauT, typename LaySourceT, typename LevIncT,
          typename LevDecT, typename SfcEmisT, typename SfcSrcT, typename RadnUpT, typename RadnDnT>
void lw_solver_noscat(int ncol, int nlay, int ngpt, bool top_at_1, DT const &D, WeightsT const &weights, int weight_ind, TauT const &tau,
                      LaySourceT const &lay_source, LevIncT const &lev_source_inc, LevDecT const &lev_source_dec,
                      SfcEmisT const &sfc_emis, SfcSrcT const &sfc_src, RadnUpT const &radn_up, RadnDnT const &radn_dn) {
  using RealT   = typename TauT::non_const_value_type;
  using LayoutT = typename TauT::array_layout;
  using DeviceT = typename TauT::device_type;
  using pool = conv::MemPoolSingleton<RealT, DeviceT>;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  using ureal2d_t = conv::Unmanaged<Kokkos::View<RealT**, LayoutT, DeviceT>>;
  using real3d_t  = Kokkos::View<RealT***, LayoutT, DeviceT>;
  using ureal3d_t = conv::Unmanaged<real3d_t>;

  const int dsize1 = ncol*nlay*ngpt;
  const int dsize2 = ncol*ngpt;
  RealT* data = pool::template alloc_raw<RealT>(dsize1*4 + dsize2*2), *dcurr=data;
  ureal3d_t tau_loc   (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t trans     (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t source_dn (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal3d_t source_up (dcurr,ncol,nlay,ngpt); dcurr += dsize1;
  ureal2d_t source_sfc(dcurr,ncol,     ngpt); dcurr += dsize2;
  ureal2d_t sfc_albedo(dcurr,ncol,     ngpt); dcurr += dsize2;

  RealT tau_thresh = sqrt( std::numeric_limits<RealT>::epsilon() );

  RealT constexpr pi = M_PI;

  // Which way is up?
  // Level Planck sources for upward and downward radiation
  // When top_at_1, lev_source_up => lev_source_dec
  //                lev_source_dn => lev_source_inc, and vice-versa
  int top_level;
  real3d_t lev_source_up;
  real3d_t lev_source_dn;
  if (top_at_1) {
    top_level = 0;
    // Recall below that equating two arrays is like assigning pointers in Fortran. No data is copied.
    // The LHS just uses the same data pointer as the RHS so that changing one's data changes the other's as well.
    lev_source_up = lev_source_dec;
    lev_source_dn = lev_source_inc;
  } else {
    top_level = nlay;
    // Recall below that equating two arrays is like assigning pointers in Fortran. No data is copied.
    // The LHS just uses the same data pointer as the RHS so that changing one's data changes the other's as well.
    lev_source_up = lev_source_inc;
    lev_source_dn = lev_source_dec;
  }

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
    // Transport is for intensity
    //   convert flux at top of domain to intensity assuming azimuthal isotropy
    radn_dn(icol,top_level,igpt) = radn_dn(icol,top_level,igpt)/(2. * pi * weights(weight_ind));

    // Surface albedo, surface source function
    sfc_albedo(icol,igpt) = 1. - sfc_emis(icol,igpt);
    source_sfc(icol,igpt) = sfc_emis(icol,igpt) * sfc_src(icol,igpt);
  }));

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  Kokkos::Array<int, 3> dims3_ngpt_nlay_ncol = {ncol,nlay,ngpt};
  const int dims3_ngpt_nlay_ncol_tot = ngpt * nlay * ncol;
  //TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ngpt,nlay,ncol}) , KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
  TIMED_KERNEL(Kokkos::parallel_for( dims3_ngpt_nlay_ncol_tot , KOKKOS_LAMBDA (int idx) {
    int icol, ilay, igpt;
    conv::unflatten_idx_left(idx, dims3_ngpt_nlay_ncol, icol, ilay, igpt);
    // Optical path and transmission, used in source function and transport calculations
    tau_loc(icol,ilay,igpt) = tau(icol,ilay,igpt)*D(icol,igpt);
    trans  (icol,ilay,igpt) = exp(-tau_loc(icol,ilay,igpt));

    lw_source_noscat_stencil(ncol, nlay, ngpt, icol, ilay, igpt,
                             lay_source, lev_source_up, lev_source_dn,
                             tau_loc, trans,
                             source_dn, source_up, tau_thresh);
  }));

  // Transport
  lw_transport_noscat(ncol, nlay, ngpt, top_at_1,
                      tau_loc, trans, sfc_albedo, source_dn, source_up, source_sfc,
                      radn_up, radn_dn);

  // Convert intensity to flux assuming azimuthal isotropy and quadrature weight
  // do igpt = 1, ngpt
  //   do ilev = 1, nlay+1
  //     do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template getrl<3>({ncol,nlay+1,ngpt}) , KOKKOS_LAMBDA (int icol, int ilev, int igpt) {
    radn_dn(icol,ilev,igpt) = 2. * pi * weights(weight_ind) * radn_dn(icol,ilev,igpt);
    radn_up(icol,ilev,igpt) = 2. * pi * weights(weight_ind) * radn_up(icol,ilev,igpt);
  }));

  pool::dealloc(data, dcurr - data);
}

// LW transport, no scattering, multi-angle quadrature
//   Users provide a set of weights and quadrature angles
//   Routine sums over single-angle solutions for each sets of angles/weights
template <typename DsT, typename WeightsT, typename TauT, typename LaySourceT, typename LevIncT,
          typename LevDecT, typename SfcEmisT, typename SfcSrcT, typename FluxUpT, typename FluxDnT>
void lw_solver_noscat_GaussQuad(int ncol, int nlay, int ngpt, bool top_at_1, int nmus, DsT const &Ds, WeightsT const &weights,
                                TauT const &tau, LaySourceT const &lay_source, LevIncT const &lev_source_inc, LevDecT const &lev_source_dec,
                                SfcEmisT const &sfc_emis, SfcSrcT const &sfc_src, FluxUpT const &flux_up, FluxDnT const &flux_dn) {
  // Local variables
  using RealT   = typename TauT::non_const_value_type;
  using LayoutT = typename TauT::array_layout;
  using DeviceT = typename TauT::device_type;
  using pool = conv::MemPoolSingleton<RealT, DeviceT>;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  using ureal2d_t = conv::Unmanaged<Kokkos::View<RealT**, LayoutT, DeviceT>>;
  using ureal3d_t = conv::Unmanaged<Kokkos::View<RealT***, LayoutT, DeviceT>>;

  const int dsize1 = ncol*(nlay+1)*ngpt;
  const int dsize2 = ncol*ngpt;
  RealT* data = pool::template alloc_raw<RealT>(dsize1*2 + dsize2*2), *dcurr=data;
  ureal3d_t radn_dn (dcurr,ncol,nlay+1,ngpt); dcurr += dsize1;
  ureal3d_t radn_up (dcurr,ncol,nlay+1,ngpt); dcurr += dsize1;
  ureal2d_t Ds_ncol (dcurr,ncol,       ngpt); dcurr += dsize2;
  ureal2d_t flux_top(dcurr,ncol,       ngpt); dcurr += dsize2;

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
    Ds_ncol(icol, igpt) = Ds(0);
  }));

  lw_solver_noscat(ncol, nlay, ngpt,
                   top_at_1, Ds_ncol, weights, 0, tau,
                   lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src,
                   flux_up, flux_dn);
  //
  // For more than one angle use local arrays
  int top_level = conv::merge(0, nlay, top_at_1);

  // do igpt = 1, ngpt
  //   do icol = 1, ncol
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
    flux_top(icol,igpt) = flux_dn(icol,top_level,igpt);
  }));

  apply_BC(ncol, nlay, ngpt, top_at_1, flux_top, radn_dn);

  for (int imu=1; imu<nmus; imu++) {
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      Ds_ncol(icol, igpt) = Ds(imu);
    }));

    lw_solver_noscat(ncol, nlay, ngpt,
                     top_at_1, Ds_ncol, weights, imu, tau,
                     lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src,
                     radn_up, radn_dn);

    // do igpt = 1, ngpt
    //   do ilev = 1, nlay+1
    //     do icol = 1, ncol
    TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ngpt,nlay+1,ncol}) , KOKKOS_LAMBDA (int igpt, int ilev, int icol) {
      flux_up(icol,ilev,ngpt) += radn_up(icol,ilev,ngpt);
      flux_dn(icol,ilev,ngpt) += radn_dn(icol,ilev,ngpt);
    }));

  } // imu

  pool::dealloc(data, dcurr - data);
}

#endif
