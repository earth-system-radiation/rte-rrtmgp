
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
//    atmospheric optical properties on a spectral grid
//    information about vertical ordering
//    boundary conditions
//      solar zenith angle, spectrally-resolved incident colimated flux, surface albedos for direct and diffuse radiation
//    optionally, a boundary condition for incident diffuse radiation
//
// It is the user's responsibility to ensure that boundary conditions (incident fluxes, surface albedos) are on the same
//   spectral grid as the optical properties.
//
// Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
//   whatever summary the user needs.
//
// The routine does error checking and choses which lower-level kernel to invoke based on
//   what kinds of optical properties are supplied
//
// -------------------------------------------------------------------------------------------------

#ifdef RRTMGP_ENABLE_YAKL
template <class FluxesType>
void rte_sw(OpticalProps2str const &atmos, bool top_at_1, real1d const &mu0, real2d const &inc_flux,
            real2d const &sfc_alb_dir, real2d const &sfc_alb_dif, FluxesType &fluxes, real2d const &inc_flux_dif=real2d());

template <class FluxesType>
void rte_sw(OpticalProps2str const &atmos, bool top_at_1, real1d const &mu0, real2d const &inc_flux,
            real2d const &sfc_alb_dir, real2d const &sfc_alb_dif, FluxesType &fluxes, real2d const &inc_flux_dif) {
  using yakl::intrinsics::size;
  using yakl::intrinsics::allocated;
  using yakl::intrinsics::any;
  using yakl::componentwise::operator<;
  using yakl::componentwise::operator>;

  real3d gpt_flux_up;
  real3d gpt_flux_dn;
  real3d gpt_flux_dir;
  real2d sfc_alb_dir_gpt;
  real2d sfc_alb_dif_gpt;
  int ncol  = atmos.get_ncol();
  int nlay  = atmos.get_nlay();
  int ngpt  = atmos.get_ngpt();
  int nband = atmos.get_nband();

  // Error checking -- consistency of sizes and validity of values
  if (! fluxes.are_desired()) { stoprun("rte_sw: no space allocated for fluxes"); }

  if (size(mu0,1) != ncol) { stoprun("rte_sw: mu0 inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(mu0 < 0.) || any(mu0 > 1.)) { stoprun("rte_sw: one or more mu0 <= 0 or > 1"); }
  #endif

  if (size(inc_flux,1) != ncol || size(inc_flux,2) != ngpt) { stoprun("rte_sw: inc_flux inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(inc_flux < 0.)) { stoprun("rte_sw: one or more inc_flux < 0"); }
  #endif
  if (allocated(inc_flux_dif)) {
    if (size(inc_flux_dif,1) != ncol || size(inc_flux_dif,2) != ngpt) { stoprun("rte_sw: inc_flux_dif inconsistently sized"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(inc_flux_dif < 0.)) { stoprun("rte_sw: one or more inc_flux_dif < 0"); }
    #endif
  }

  if (size(sfc_alb_dir,1) != nband || size(sfc_alb_dir,2) != ncol) { stoprun("rte_sw: sfc_alb_dir inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(sfc_alb_dir < 0.) || any(sfc_alb_dir > 1.)) { stoprun("rte_sw: sfc_alb_dir out of bounds [0,1]"); }
  #endif
  if (size(sfc_alb_dif,1) != nband || size(sfc_alb_dif,2) != ncol) { stoprun("rte_sw: sfc_alb_dif inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(sfc_alb_dif < 0.) || any(sfc_alb_dif > 1.)) { stoprun("rte_sw: sfc_alb_dif out of bounds [0,1]"); }
  #endif

  gpt_flux_up  = real3d("gpt_flux_up" ,ncol, nlay+1, ngpt);
  gpt_flux_dn  = real3d("gpt_flux_dn" ,ncol, nlay+1, ngpt);
  gpt_flux_dir = real3d("gpt_flux_dir",ncol, nlay+1, ngpt);
  sfc_alb_dir_gpt = real2d("sfc_alb_dir_gpt",ncol, ngpt);
  sfc_alb_dif_gpt = real2d("sfc_alb_dif_gpt",ncol, ngpt);
  // Lower boundary condition -- expand surface albedos by band to gpoints
  //   and switch dimension ordering
  expand_and_transpose(atmos, sfc_alb_dir, sfc_alb_dir_gpt);
  expand_and_transpose(atmos, sfc_alb_dif, sfc_alb_dif_gpt);

  // Compute the radiative transfer...
  // Apply boundary conditions
  //   On input flux_dn is the diffuse component; the last action in each solver is to add
  //   direct and diffuse to represent the total, consistent with the LW
  apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, mu0, gpt_flux_dir);
  if (allocated(inc_flux_dif)) {
    apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif,  gpt_flux_dn );
  } else {
    apply_BC(ncol, nlay, ngpt, top_at_1,                gpt_flux_dn );
  }

  // two-stream calculation with scattering
  atmos.validate();
  sw_solver_2stream(ncol, nlay, ngpt, top_at_1,
                    atmos.tau, atmos.ssa, atmos.g, mu0,
                    sfc_alb_dir_gpt, sfc_alb_dif_gpt,
                    gpt_flux_up, gpt_flux_dn, gpt_flux_dir);

  // ...and reduce spectral fluxes to desired output quantities
  fluxes.reduce(gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir);
}
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
template <typename RealT, typename LayoutT, typename DeviceT, typename Mu0T, typename IncFluxT,
          typename SfcDirT, typename SfcDifT, class FluxesType,
          typename IncFluxDifT = typename Kokkos::View<RealT**, LayoutT, DeviceT> >
void rte_sw(OpticalProps2strK<RealT, LayoutT, DeviceT> const &atmos, bool top_at_1,
            Mu0T const &mu0, IncFluxT const &inc_flux,
            SfcDirT const &sfc_alb_dir, SfcDifT const &sfc_alb_dif,
            FluxesType &fluxes, IncFluxDifT const &inc_flux_dif=IncFluxDifT())
{
  using pool = conv::MemPoolSingleton<RealT, LayoutT, DeviceT>;
  using ureal2d_t = conv::Unmanaged<Kokkos::View<RealT**,  LayoutT, DeviceT>>;
  using ureal3d_t = conv::Unmanaged<Kokkos::View<RealT***, LayoutT, DeviceT>>;

  const int ncol  = atmos.get_ncol();
  const int nlay  = atmos.get_nlay();
  const int ngpt  = atmos.get_ngpt();
  const int nband = atmos.get_nband();

  const int dsize1 = ncol * (nlay+1) * ngpt;
  const int dsize2 = ncol * ngpt;
  RealT* data = pool::template alloc_raw<RealT>(dsize1*3 + dsize2*2), *dcurr = data;
  ureal3d_t gpt_flux_up    (dcurr,ncol, nlay+1, ngpt); dcurr += dsize1;
  ureal3d_t gpt_flux_dn    (dcurr,ncol, nlay+1, ngpt); dcurr += dsize1;
  ureal3d_t gpt_flux_dir   (dcurr,ncol, nlay+1, ngpt); dcurr += dsize1;
  ureal2d_t sfc_alb_dir_gpt(dcurr,ncol, ngpt);         dcurr += dsize2;
  ureal2d_t sfc_alb_dif_gpt(dcurr,ncol, ngpt);         dcurr += dsize2;

  // Error checking -- consistency of sizes and validity of values
  if (! fluxes.are_desired()) { stoprun("rte_sw: no space allocated for fluxes"); }

  if (mu0.extent(0) != ncol) { stoprun("rte_sw: mu0 inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(mu0 < 0.) || any(mu0 > 1.)) { stoprun("rte_sw: one or more mu0 <= 0 or > 1"); }
  #endif

  if (inc_flux.extent(0) != ncol || inc_flux.extent(1) != ngpt) { stoprun("rte_sw: inc_flux inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(inc_flux < 0.)) { stoprun("rte_sw: one or more inc_flux < 0"); }
  #endif
  if (inc_flux_dif.is_allocated()) {
    if (inc_flux_dif.extent(0) != ncol || inc_flux_dif.extent(1) != ngpt) { stoprun("rte_sw: inc_flux_dif inconsistently sized"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(inc_flux_dif < 0.)) { stoprun("rte_sw: one or more inc_flux_dif < 0"); }
    #endif
  }

  if (sfc_alb_dir.extent(0) != nband || sfc_alb_dir.extent(1) != ncol) { stoprun("rte_sw: sfc_alb_dir inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(sfc_alb_dir < 0.) || any(sfc_alb_dir > 1.)) { stoprun("rte_sw: sfc_alb_dir out of bounds [0,1]"); }
  #endif
  if (sfc_alb_dif.extent(0) != nband || sfc_alb_dif.extent(1) != ncol) { stoprun("rte_sw: sfc_alb_dif inconsistently sized"); }
  #ifdef RRTMGP_EXPENSIVE_CHECKS
    if (any(sfc_alb_dif < 0.) || any(sfc_alb_dif > 1.)) { stoprun("rte_sw: sfc_alb_dif out of bounds [0,1]"); }
  #endif

  // Lower boundary condition -- expand surface albedos by band to gpoints
  //   and switch dimension ordering
  expand_and_transpose(atmos, sfc_alb_dir, sfc_alb_dir_gpt);
  expand_and_transpose(atmos, sfc_alb_dif, sfc_alb_dif_gpt);

  // Compute the radiative transfer...
  // Apply boundary conditions
  //   On input flux_dn is the diffuse component; the last action in each solver is to add
  //   direct and diffuse to represent the total, consistent with the LW
  apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux, mu0, gpt_flux_dir);
  if (inc_flux_dif.is_allocated()) {
    apply_BC(ncol, nlay, ngpt, top_at_1, inc_flux_dif,  gpt_flux_dn );
  } else {
    apply_BC(ncol, nlay, ngpt, top_at_1,                gpt_flux_dn );
  }

  // two-stream calculation with scattering
  atmos.validate();
  sw_solver_2stream(ncol, nlay, ngpt, top_at_1,
                    atmos.tau, atmos.ssa, atmos.g, mu0,
                    sfc_alb_dir_gpt, sfc_alb_dif_gpt,
                    gpt_flux_up, gpt_flux_dn, gpt_flux_dir);

  // ...and reduce spectral fluxes to desired output quantities
  fluxes.reduce(gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir);

  pool::dealloc(data, dcurr - data);
}
#endif
