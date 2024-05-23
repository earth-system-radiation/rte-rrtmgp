
#pragma once

#include "rrtmgp_const.h"
#include "rrtmgp_conversion.h"

#ifdef RRTMGP_ENABLE_YAKL
// increment 2stream by 2stream
void inc_2stream_by_2stream_bybnd(int ncol, int nlay, int ngpt,
                                  real3d const &tau1, real3d const &ssa1, real3d const &g1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2,
                                  int nbnd, int2d const &gpt_lims);



// Incrementing when the second set of optical properties is defined at lower spectral resolution
//   (e.g. by band instead of by gpoint)
void inc_1scalar_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &tau2,
                                  int nbnd, int2d const &gpt_lims);



// Delta-scale
//   f = g*g
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, real3d const &tau, real3d const &ssa, real3d const &g);


// Delta-scaling, provided only for two-stream properties at present
// -------------------------------------------------------------------------------------------------
// Delta-scale
//   user-provided value of f (forward scattering)
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, real3d const &tau, real3d const &ssa, real3d const &g, real3d const &f);



// Addition of optical properties: the first set are incremented by the second set.
//
//   There are three possible representations of optical properties (scalar = optical depth only;
//   two-stream = tau, single-scattering albedo, and asymmetry factor g, and
//   n-stream = tau, ssa, and phase function moments p.) Thus we need nine routines, three for
//   each choice of representation on the left hand side times three representations of the
//   optical properties to be added.
//
//   There are two sets of these nine routines. In the first the two sets of optical
//   properties are defined at the same spectral resolution. There is also a set of routines
//   to add properties defined at lower spectral resolution to a set defined at higher spectral
//   resolution (adding properties defined by band to those defined by g-point)
void increment_1scalar_by_1scalar(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &tau2);



// increment 1scalar by 2stream
void increment_1scalar_by_2stream(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &tau2, real3d const &ssa2);



// increment 1scalar by nstream
void increment_1scalar_by_nstream(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &tau2, real3d const &ssa2);



// increment 2stream by 1scalar
void increment_2stream_by_1scalar(int ncol, int nlay, int ngpt, real3d const &tau1, real3d &ssa1, real3d const &tau2);



// increment 2stream by 2stream
void increment_2stream_by_2stream(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &ssa1, real3d const &g1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2);


// increment 2stream by nstream
void increment_2stream_by_nstream(int ncol, int nlay, int ngpt, int nmom2, real3d const &tau1, real3d const &ssa1,
                                  real3d &g1, real3d const &tau2, real3d const &ssa2, real4d const &p2);



// increment nstream by 1scalar
void increment_nstream_by_1scalar(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &ssa1, real3d const &tau2);



// increment nstream by 2stream
void increment_nstream_by_2stream(int ncol, int nlay, int ngpt, int nmom1, real3d const &tau1, real3d const &ssa1,
                                  real4d const &p1, real3d const &tau2, real3d const &ssa2, real3d const &g2);



// increment nstream by nstream
void increment_nstream_by_nstream(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real3d const &tau1,
                                  real3d const &ssa1, real4d const &p1, real3d const &tau2, real3d const &ssa2, real4d const &p2);


// increment 1scalar by 2stream
void inc_1scalar_by_2stream_bybnd(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &tau2,
                                  real3d const &ssa2, int nbnd, int2d const &gpt_lims);



// increment 1scalar by nstream
void inc_1scalar_by_nstream_bybnd(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &tau2,
                                  real3d const &ssa2, int nbnd, int2d const &gpt_lims);



// increment 2stream by 1scalar
void inc_2stream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &ssa1,
                                  real3d const &tau2, int nbnd, int2d const &gpt_lims);


// increment 2stream by nstream
void inc_2stream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom2, real3d const &tau1, real3d const &ssa1,
                                  real3d const &g1, real3d const &tau2, real3d const &ssa2, real4d const &p2,
                                  int nbnd, int2d const &gpt_lims);



// increment nstream by 1scalar
void inc_nstream_by_1scalar_bybnd(int ncol, int nlay, int ngpt, real3d const &tau1, real3d const &ssa1, real3d const &tau2,
                                  int nbnd, int2d const &gpt_lims);



// increment nstream by 2stream
void inc_nstream_by_2stream_bybnd(int ncol, int nlay, int ngpt, int nmom1, real3d const &tau1, real3d const &ssa1, real4d const &p1,
                                  real3d const &tau2, real3d const &ssa2, real3d const &g2, int nbnd,
                                  int2d const &gpt_lims);



// increment nstream by nstream
void inc_nstream_by_nstream_bybnd(int ncol, int nlay, int ngpt, int nmom1, int nmom2, real3d const &tau1, real3d const &ssa1, real4d const &p1,
                                  real3d const &tau2, real3d const &ssa2, real4d const &p2, int nbnd, int2d const &gpt_lims);



// Subsetting, meaning extracting some portion of the 3D domain
void extract_subset_dim1_3d(int ncol, int nlay, int ngpt, real3d const &array_in, int colS, int colE, real3d const &array_out);



void extract_subset_dim2_4d(int nmom, int ncol, int nlay, int ngpt, real3d const &array_in, int colS, int colE, real3d const &array_out);



// Extract the absorption optical thickness which requires mulitplying by 1 - ssa
void extract_subset_absorption_tau(int ncol, int nlay, int ngpt, real3d const &tau_in, real3d const &ssa_in, int colS, int colE,
                                   real3d const &tau_out);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS

// increment 2stream by 2stream
template <typename Tau1T, typename Ssa1T, typename G1T,
          typename Tau2T, typename Ssa2T, typename G2T,
          typename GptLimsT>
void inc_2stream_by_2stream_bybnd(int ncol, int nlay, int ngpt,
                                  Tau1T const &tau1, Ssa1T const &ssa1, G1T const &g1,
                                  Tau2T const &tau2, Ssa2T const &ssa2, G2T const &g2,
                                  int nbnd, GptLimsT const &gpt_lims) {
  using RealT = typename Tau1T::non_const_value_type;
  constexpr RealT eps = 3*std::numeric_limits<RealT>::min();

  // do igpt = 1 , ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  Kokkos::parallel_for( conv::get_mdrp<4>({nbnd,ngpt,nlay,ncol}) , KOKKOS_LAMBDA (int ibnd, int igpt, int ilay, int icol) {
    if (igpt >= gpt_lims(0,ibnd) && igpt <= gpt_lims(1,ibnd) ) {
      // t=tau1 + tau2
      RealT tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,ibnd);
      // w=(tau1*ssa1 + tau2*ssa2) / t
      RealT tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) +
                       tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd);
      g1(icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * g1(icol,ilay,igpt) +
                            tau2(icol,ilay,ibnd) * ssa2(icol,ilay,ibnd) * g2(icol,ilay,ibnd)) / Kokkos::fmax(eps,tauscat12);
      ssa1(icol,ilay,igpt) = tauscat12 / Kokkos::fmax(eps,tau12);
      tau1(icol,ilay,igpt) = tau12;
    }
  });
}

// Addition of optical properties: the first set are incremented by the second set.
//
//   There are three possible representations of optical properties (scalar = optical depth only;
//   two-stream = tau, single-scattering albedo, and asymmetry factor g, and
//   n-stream = tau, ssa, and phase function moments p.) Thus we need nine routines, three for
//   each choice of representation on the left hand side times three representations of the
//   optical properties to be added.
//
//   There are two sets of these nine routines. In the first the two sets of optical
//   properties are defined at the same spectral resolution. There is also a set of routines
//   to add properties defined at lower spectral resolution to a set defined at higher spectral
//   resolution (adding properties defined by band to those defined by g-point)
template <typename Tau1T, typename Tau2T>
void increment_1scalar_by_1scalar(int ncol, int nlay, int ngpt, Tau1T const &tau1, Tau2T const &tau2) {
  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  Kokkos::parallel_for( conv::get_mdrp<3>({ngpt,nlay,ncol}), KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    tau1(icol,ilay,igpt) = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
  });
  //std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}

// Incrementing when the second set of optical properties is defined at lower spectral resolution
//   (e.g. by band instead of by gpoint)
template <typename Tau1T, typename Tau2T, typename GptLimsT>
void inc_1scalar_by_1scalar_bybnd(int ncol, int nlay, int ngpt, Tau1T const &tau1, Tau2T const &tau2,
                                  int nbnd, GptLimsT const &gpt_lims) {
  // do igpt = 1 , ngpt
  //   do ilay = 1 , nlay
  //     do icol = 1 , ncol
  Kokkos::parallel_for( conv::get_mdrp<3>({ngpt,nlay,ncol}) , KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    for (int ibnd=0; ibnd<nbnd; ibnd++) {
      if (igpt >= gpt_lims(0,ibnd) && igpt <= gpt_lims(1,ibnd) ) {
        tau1(icol,ilay,igpt) += tau2(icol,ilay,ibnd);
      }
    }
  });
}

// Delta-scale
//   f = g*g
template <typename TauT, typename SsaT, typename GT>
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, TauT const &tau, SsaT const &ssa, GT const &g) {
  using RealT = typename TauT::non_const_value_type;

  constexpr RealT eps = 3*std::numeric_limits<RealT>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  Kokkos::parallel_for( conv::get_mdrp<3>({ngpt,nlay,ncol}) , KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    if (tau(icol,ilay,igpt) > eps) {
      RealT f  = g  (icol,ilay,igpt) * g  (icol,ilay,igpt);
      RealT wf = ssa(icol,ilay,igpt) * f;
      tau(icol,ilay,igpt) = (1. - wf) * tau(icol,ilay,igpt);
      ssa(icol,ilay,igpt) = (ssa(icol,ilay,igpt) - wf) / (1.0 - wf);
      g  (icol,ilay,igpt) = (g  (icol,ilay,igpt) -  f) / (1.0 -  f);
    }
  });
}


// Delta-scaling, provided only for two-stream properties at present
// -------------------------------------------------------------------------------------------------
// Delta-scale
//   user-provided value of f (forward scattering)
template <typename TauT, typename SsaT, typename GT, typename FT>
void delta_scale_2str_kernel(int ncol, int nlay, int ngpt, TauT const &tau, SsaT const &ssa, GT const &g, FT const &f) {
  using RealT = typename TauT::non_const_value_type;
  constexpr RealT eps = 3*std::numeric_limits<RealT>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  Kokkos::parallel_for( conv::get_mdrp<3>({ngpt,nlay,ncol}) , KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    if (tau(icol,ilay,igpt) > eps) {
      RealT wf = ssa(icol,ilay,igpt) * f(icol,ilay,igpt);
      tau(icol,ilay,igpt) = (1. - wf) * tau(icol,ilay,igpt);
      ssa(icol,ilay,igpt) = (ssa(icol,ilay,igpt) - wf) /  (1.0 - wf);
      g  (icol,ilay,igpt) = (g  (icol,ilay,igpt) - f(icol,ilay,igpt)) / (1. - f(icol,ilay,igpt));
    }
  });
  std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}

// increment 2stream by 2stream
template <typename Tau1T, typename Ssa1T, typename G1T,
          typename Tau2T, typename Ssa2T, typename G2T>
void increment_2stream_by_2stream(int ncol, int nlay, int ngpt, Tau1T const &tau1, Ssa1T const &ssa1, G1T const &g1,
                                  Tau2T const &tau2, Ssa2T const &ssa2, G2T const &g2) {
  using RealT = typename Tau1T::non_const_value_type;
  constexpr RealT eps = 3*std::numeric_limits<RealT>::min();

  // do igpt = 1, ngpt
  //   do ilay = 1, nlay
  //     do icol = 1, ncol
  Kokkos::parallel_for( conv::get_mdrp<3>({ngpt,nlay,ncol}) , KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    // t=tau1 + tau2
    RealT tau12 = tau1(icol,ilay,igpt) + tau2(icol,ilay,igpt);
    // w=(tau1*ssa1 + tau2*ssa2) / t
    RealT tauscat12 = tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) +
                     tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt);
    g1(icol,ilay,igpt) = (tau1(icol,ilay,igpt) * ssa1(icol,ilay,igpt) * g1(icol,ilay,igpt) +
                          tau2(icol,ilay,igpt) * ssa2(icol,ilay,igpt) * g2(icol,ilay,igpt))
                           / Kokkos::fmax(eps,tauscat12);
    ssa1(icol,ilay,igpt) = tauscat12 / Kokkos::fmax(eps,tau12);
    tau1(icol,ilay,igpt) = tau12;
  });
  //std::cout << "WARNING: THIS ISN'T TESTED: " << __FILE__ << ": " << __LINE__ << "\n";
}

#endif
