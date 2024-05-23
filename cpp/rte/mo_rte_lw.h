
#pragma once

#include "rrtmgp_const.h"
#include "mo_optical_props.h"
#include "mo_source_functions.h"
#include "expand_and_transpose.h"
#include "mo_rte_solver_kernels.h"

// This code is part of Radiative Transfer for Energetics (RTE)
//
// Contacts: Robert Pincus and Eli Mlawer
// email:  rrtmgp@aer.com
//
// Copyright 2015-2018,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
// -------------------------------------------------------------------------------------------------
//
//  Contains a single routine to compute direct and diffuse fluxes of solar radiation given
//    atmospheric optical properties, spectrally-resolved
//    information about vertical ordering
//    internal Planck source functions, defined per g-point on the same spectral grid at the atmosphere
//    boundary conditions: surface emissivity defined per band
//    optionally, a boundary condition for incident diffuse radiation
//    optionally, an integer number of angles at which to do Gaussian quadrature if scattering is neglected
//
// If optical properties are supplied via class ty_optical_props_1scl (absorption optical thickenss only)
//    then an emission/absorption solver is called
//    If optical properties are supplied via class ty_optical_props_2str fluxes are computed via
//    two-stream calculations and adding.
//
// It is the user's responsibility to ensure that emissivity is on the same
//   spectral grid as the optical properties.
//
// Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
//   whatever summary the user needs.
//
// The routine does error checking and choses which lower-level kernel to invoke based on
//   what kinds of optical properties are supplied
//
// -------------------------------------------------------------------------------------------------


// Interface using only optical properties and source functions as inputs; fluxes as outputs.
#ifdef RRTMGP_ENABLE_YAKL
template <class FluxesType>
void rte_lw(int max_gauss_pts, real2d const &gauss_Ds, real2d const &gauss_wts, OpticalProps1scl const &optical_props,
            bool top_at_1, SourceFuncLW const &sources, real2d const &sfc_emis,
            FluxesType &fluxes, real2d const &inc_flux=real2d(), int n_gauss_angles=-1);

template <class FluxesType>
void rte_lw(int max_gauss_pts, real2d const &gauss_Ds, real2d const &gauss_wts, OpticalProps1scl const &optical_props,
            bool top_at_1, SourceFuncLW const &sources, real2d const &sfc_emis,
            FluxesType &fluxes, real2d const &inc_flux, int n_gauss_angles) {
  using yakl::intrinsics::size;
  using yakl::intrinsics::allocated;
  using yakl::intrinsics::any;
  using yakl::componentwise::operator<;
  using yakl::componentwise::operator>;
  using yakl::fortran::parallel_for;
  using yakl::fortran::SimpleBounds;

  real3d gpt_flux_up;
  real3d gpt_flux_dn;
  real2d sfc_emis_gpt;

  // Error checking
  //   if inc_flux is present it has the right dimensions, is positive definite
  int ncol  = optical_props.get_ncol();
  int nlay  = optical_props.get_nlay();
  int ngpt  = optical_props.get_ngpt();
  int nband = optical_props.get_nband();

  // Error checking -- consistency of sizes and validity of values
  if (! fluxes.are_desired()) { stoprun("rte_lw: no space allocated for fluxes"); }

  // Source functions
  if (sources.get_ncol() != ncol || sources.get_nlay() != nlay || sources.get_ngpt() != ngpt) {
    stoprun("rte_lw: sources and optical properties inconsistently sized");
  }

  // Surface emissivity
  if (size(sfc_emis,1) != nband || size(sfc_emis,2) != ncol) { stoprun("rte_lw: sfc_emis inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(sfc_emis < 0.) || any(sfc_emis > 1.)) { stoprun("rte_lw: sfc_emis has values < 0 or > 1"); }
  #endif

  // Incident flux, if present
  if (allocated(inc_flux)) {
    if (size(inc_flux,1) != ncol | size(inc_flux,2) != ngpt) { stoprun("rte_lw: inc_flux inconsistently sized"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(inc_flux < 0.)) { stoprun("rte_lw: inc_flux has values < 0"); }
    #endif
  }

  // Number of quadrature points for no-scattering calculation
  int n_quad_angs = 1;
  if ( n_gauss_angles != -1 ) {
    if (n_gauss_angles > max_gauss_pts) {
      stoprun("rte_lw: asking for too many quadrature points for no-scattering calculation");
    }
    if (n_gauss_angles < 1) {
      stoprun("rte_lw: have to ask for at least one quadrature point for no-scattering calculation");
    }
    n_quad_angs = n_gauss_angles;
  }

  // Ensure values of tau, ssa, and g are reasonable
  optical_props.validate();

  // Lower boundary condition -- expand surface emissivity by band to gpoints
  gpt_flux_up  = real3d("gpt_flux_up" ,ncol,nlay+1,ngpt);
  gpt_flux_dn  = real3d("gpt_flux_dn" ,ncol,nlay+1,ngpt);
  sfc_emis_gpt = real2d("sfc_emis_gpt",ncol       ,ngpt);
  expand_and_transpose(optical_props, sfc_emis, sfc_emis_gpt);

  //   Upper boundary condition
  if (allocated(inc_flux)) {
    apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, gpt_flux_dn);
  } else {
    // Default is zero incident diffuse flux
    apply_BC(ncol, nlay, ngpt, top_at_1          , gpt_flux_dn);
  }

  // Compute the radiative transfer...
  // No scattering two-stream calculation
  optical_props.validate();
  real1d tmp_Ds ("tmp_Ds" ,n_quad_angs);
  real1d tmp_wts("tmp_wts",n_quad_angs);
  // for (int i=1 ; i <= n_quad_angs ; i++) {
  parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<1>(n_quad_angs) , YAKL_LAMBDA (int i) {
    tmp_Ds (i) = gauss_Ds (i,n_quad_angs);
    tmp_wts(i) = gauss_wts(i,n_quad_angs);
  });
  lw_solver_noscat_GaussQuad(ncol, nlay, ngpt, top_at_1, n_quad_angs, tmp_Ds, tmp_wts, optical_props.tau,
                             sources.lay_source, sources.lev_source_inc, sources.lev_source_dec,
                             sfc_emis_gpt, sources.sfc_source, gpt_flux_up, gpt_flux_dn);
  // ...and reduce spectral fluxes to desired output quantities
  fluxes.reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1);
}
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
// Interface using only optical properties and source functions as inputs; fluxes as outputs.
template <typename GaussDsT, typename GaussWtsT, typename RealT, typename LayoutT, typename DeviceT, typename SfcEmisT,
          class FluxesType, typename IncFluxT=typename Kokkos::View<RealT**, LayoutT, DeviceT> >
void rte_lw(int max_gauss_pts, GaussDsT const &gauss_Ds, GaussWtsT const &gauss_wts,
            OpticalProps1sclK<RealT, LayoutT, DeviceT> const &optical_props,
            bool top_at_1, SourceFuncLWK<RealT, LayoutT, DeviceT> const &sources, SfcEmisT const &sfc_emis,
            FluxesType &fluxes, IncFluxT const &inc_flux=IncFluxT(), int n_gauss_angles=-1)
{
  using pool = conv::MemPoolSingleton;
  using real1d_t = Kokkos::View<RealT*,   LayoutT, DeviceT>;
  using real2d_t = Kokkos::View<RealT**,  LayoutT, DeviceT>;
  using real3d_t = Kokkos::View<RealT***, LayoutT, DeviceT>;

  const int ncol  = optical_props.get_ncol();
  const int nlay  = optical_props.get_nlay();
  const int ngpt  = optical_props.get_ngpt();
  const int nband = optical_props.get_nband();

  // Number of quadrature points for no-scattering calculation
  int n_quad_angs = 1;
  if ( n_gauss_angles != -1 ) {
    if (n_gauss_angles > max_gauss_pts) {
      stoprun("rte_lw: asking for too many quadrature points for no-scattering calculation");
    }
    if (n_gauss_angles < 1) {
      stoprun("rte_lw: have to ask for at least one quadrature point for no-scattering calculation");
    }
    n_quad_angs = n_gauss_angles;
  }

  const int dsize1 = ncol * (nlay+1) * ngpt;
  const int dsize2 = ncol * ngpt;
  RealT* data = pool::alloc<real>(dsize1*2 + dsize2 + 2*n_quad_angs), *dcurr = data;
  real3d_t gpt_flux_up (dcurr,ncol,nlay+1,ngpt); dcurr += dsize1;
  real3d_t gpt_flux_dn (dcurr,ncol,nlay+1,ngpt); dcurr += dsize1;
  real2d_t sfc_emis_gpt(dcurr,ncol       ,ngpt); dcurr += dsize2;
  real1d_t tmp_Ds      (dcurr,n_quad_angs); dcurr += n_quad_angs;
  real1d_t tmp_wts     (dcurr,n_quad_angs); dcurr += n_quad_angs;

  // Error checking
  //   if inc_flux is present it has the right dimensions, is positive definite

  // Error checking -- consistency of sizes and validity of values
  if (! fluxes.are_desired()) { stoprun("rte_lw: no space allocated for fluxes"); }

  // Source functions
  if (sources.get_ncol() != ncol || sources.get_nlay() != nlay || sources.get_ngpt() != ngpt) {
    stoprun("rte_lw: sources and optical properties inconsistently sized");
  }

  // Surface emissivity
  if (sfc_emis.extent(0) != nband || sfc_emis.extent(1) != ncol) { stoprun("rte_lw: sfc_emis inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(sfc_emis < 0.) || any(sfc_emis > 1.)) { stoprun("rte_lw: sfc_emis has values < 0 or > 1"); }
  #endif

  // Incident flux, if present
  if (inc_flux.is_allocated()) {
    if (inc_flux.extent(0) != ncol | inc_flux.extent(1) != ngpt) { stoprun("rte_lw: inc_flux inconsistently sized"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(inc_flux < 0.)) { stoprun("rte_lw: inc_flux has values < 0"); }
    #endif
  }

  // Ensure values of tau, ssa, and g are reasonable
  optical_props.validate();

  // Lower boundary condition -- expand surface emissivity by band to gpoints
  expand_and_transpose(optical_props, sfc_emis, sfc_emis_gpt);

  //   Upper boundary condition
  if (inc_flux.is_allocated()) {
    apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, gpt_flux_dn);
  } else {
    // Default is zero incident diffuse flux
    apply_BC(ncol, nlay, ngpt, top_at_1          , gpt_flux_dn);
  }

  // Compute the radiative transfer...
  // No scattering two-stream calculation
  optical_props.validate();
  // for (int i=1 ; i <= n_quad_angs ; i++) {
  Kokkos::parallel_for( n_quad_angs , KOKKOS_LAMBDA (int i) {
    tmp_Ds (i) = gauss_Ds (i,n_quad_angs - 1);
    tmp_wts(i) = gauss_wts(i,n_quad_angs - 1);
  });
  lw_solver_noscat_GaussQuad(ncol, nlay, ngpt, top_at_1, n_quad_angs, tmp_Ds, tmp_wts, optical_props.tau,
                             sources.lay_source, sources.lev_source_inc, sources.lev_source_dec,
                             sfc_emis_gpt, sources.sfc_source, gpt_flux_up, gpt_flux_dn);
  // ...and reduce spectral fluxes to desired output quantities
  fluxes.reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1);

  pool::dealloc(data, dcurr - data);
}
#endif
