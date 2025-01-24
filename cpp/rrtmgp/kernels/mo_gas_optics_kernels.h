#pragma once

#include "rrtmgp_const.h"
#include "rrtmgp_conversion.h"
#include <limits>

// This code is part of
// RRTM for GCM Applications - Parallel (RRTMGP)
//
// Eli Mlawer and Robert Pincus
// Andre Wehe and Jennifer Delamere
// email:  rrtmgp@aer.com
//
// Copyright 2015,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
//
// Description: Numeric calculations for gas optics. Absorption and Rayleigh optical depths,
//   source functions.



#ifdef RRTMGP_ENABLE_YAKL
// Compute interpolation coefficients
// for calculations of major optical depths, minor optical depths, Rayleigh,
// and Planck fractions
void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp, int2d const &flavor,
                   real1d const &press_ref_log, real1d const &temp_ref, real press_ref_log_delta, real temp_ref_min,
                   real temp_ref_delta, real press_ref_trop_log, real3d const &vmr_ref, real2d const &play,
                   real2d const &tlay, real3d const &col_gas, int2d const &jtemp, real6d const &fmajor, real5d const &fminor,
                   real4d const &col_mix, bool2d const &tropo, int4d const &jeta, int2d const &jpress);


// Combine absoprtion and Rayleigh optical depths for total tau, ssa, g
void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real3d const &tau_abs, real3d const &tau_rayleigh,
                              real3d const &tau, real3d const &ssa, real3d const &g);



//  interpolation in temperature, pressure, and eta
YAKL_INLINE real interpolate3D(real1d const &scaling, real3d const &fmajor, real4d const &k, int igpt, int1d const &jeta,
                               int jtemp, int jpress, int ngpt, int neta, int npres, int ntemp) {
  // each code block is for a different reference temperature
  real ret = scaling(1) * ( fmajor(1,1,1) * k(igpt, jeta(1)  , jpress-1, jtemp  ) +
                            fmajor(2,1,1) * k(igpt, jeta(1)+1, jpress-1, jtemp  ) +
                            fmajor(1,2,1) * k(igpt, jeta(1)  , jpress  , jtemp  ) +
                            fmajor(2,2,1) * k(igpt, jeta(1)+1, jpress  , jtemp  ) ) +
             scaling(2) * ( fmajor(1,1,2) * k(igpt, jeta(2)  , jpress-1, jtemp+1) +
                            fmajor(2,1,2) * k(igpt, jeta(2)+1, jpress-1, jtemp+1) +
                            fmajor(1,2,2) * k(igpt, jeta(2)  , jpress  , jtemp+1) +
                            fmajor(2,2,2) * k(igpt, jeta(2)+1, jpress  , jtemp+1) );
  return ret;
}



//   This function returns a single value from a subset (in gpoint) of the k table
YAKL_INLINE real interpolate2D(real2d const &fminor, real3d const &k, int igpt, int1d const &jeta, int jtemp,
                               int ngpt, int neta, int ntemp) {
  return fminor(1,1) * k(igpt, jeta(1)  , jtemp  ) +
         fminor(2,1) * k(igpt, jeta(1)+1, jtemp  ) +
         fminor(1,2) * k(igpt, jeta(2)  , jtemp+1) +
         fminor(2,2) * k(igpt, jeta(2)+1, jtemp+1);
}



// ----------------------------------------------------------
//
// One dimensional interpolation -- return all values along second table dimension
//
YAKL_INLINE void interpolate1D(real val, real offset, real delta, real2d const &table,
                               real1d const &res, int tab_d1, int tab_d2) {
  real val0 = (val - offset) / delta;
  real frac = val0 - int(val0); // get fractional part
  int index = std::min(tab_d1-1, std::max(1, (int)(val0)+1)); // limit the index range
  for (int i=1; i<=tab_d2; i++) {
    res(i) = table(index,i) + frac * (table(index+1,i) - table(index,i));
  }
}



void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres, int ntemp, int nPlanckTemp,
                           real2d const &tlay, real2d const &tlev, real1d const &tsfc, int sfc_lay, real6d const &fmajor,
                           int4d const &jeta, bool2d const &tropo, int2d const &jtemp, int2d const &jpress,
                           int1d const &gpoint_bands, int2d const &band_lims_gpt, real4d const &pfracin, real temp_ref_min,
                           real totplnk_delta, real2d const &totplnk, int2d const &gpoint_flavor, real2d const &sfc_src,
                           real3d const &lay_src, real3d const &lev_src_inc, real3d const &lev_src_dec);



// compute Rayleigh scattering optical depths
void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp,
                          int2d const &gpoint_flavor, int2d const &band_lims_gpt, real4d const &krayl, int idx_h2o,
                          real2d const &col_dry, real3d const &col_gas, real5d const &fminor, int4d const &jeta,
                          bool2d const &tropo, int2d const &jtemp, real3d const &tau_rayleigh);



// compute minor species optical depths
void gas_optical_depths_minor(int max_gpt_diff, int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                              int nminor, int nminork, int idx_h2o, int idx_tropo, int2d const &gpt_flv,
                              real3d const &kminor, int2d const &minor_limits_gpt, bool1d const &minor_scales_with_density,
                              bool1d const &scale_by_complement, int1d const &idx_minor, int1d const &idx_minor_scaling,
                              int1d const &kminor_start, real2d const &play, real2d const &tlay, real3d const &col_gas,
                              real5d const &fminor, int4d const &jeta, int2d const &layer_limits, int2d const &jtemp, real3d const &tau);



// compute minor species optical depths
void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, int2d const &gpoint_flavor, int2d const &band_lims_gpt, real4d const &kmajor,
                              real4d const &col_mix, real6d const &fmajor, int4d const &jeta, bool2d const &tropo,
                              int2d const &jtemp, int2d const &jpress, real3d const &tau);



// Compute minor and major species opitcal depth from pre-computed interpolation coefficients
//   (jeta,jtemp,jpress)
void compute_tau_absorption(int max_gpt_diff_lower, int max_gpt_diff_upper, int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp, int nminorlower,
                            int nminorklower, int nminorupper, int nminorkupper, int idx_h2o, int2d const &gpoint_flavor,
                            int2d const &band_lims_gpt, real4d const &kmajor, real3d const &kminor_lower, real3d const &kminor_upper,
                            int2d const &minor_limits_gpt_lower, int2d const &minor_limits_gpt_upper, bool1d const &minor_scales_with_density_lower,
                            bool1d const &minor_scales_with_density_upper, bool1d const &scale_by_complement_lower,
                            bool1d const &scale_by_complement_upper, int1d const &idx_minor_lower, int1d const &idx_minor_upper,
                            int1d const &idx_minor_scaling_lower, int1d const &idx_minor_scaling_upper, int1d const &kminor_start_lower,
                            int1d const &kminor_start_upper, bool2d const &tropo, real4d const &col_mix, real6d const &fmajor,
                            real5d const &fminor, real2d const &play, real2d const &tlay, real3d const &col_gas, int4d const &jeta,
                            int2d const &jtemp, int2d const &jpress, real3d const &tau, bool top_at_1);



// Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
//   using Rayleigh scattering phase function
void combine_and_reorder_nstr(int ncol, int nlay, int ngpt, int nmom, real3d const &tau_abs, real3d const &tau_rayleigh, real3d const &tau, real3d const &ssa, real4d const &p);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
//   This function returns a single value from a subset (in gpoint) of the k table
template <typename FminorT, typename KT, typename JetaT>
KOKKOS_INLINE_FUNCTION
auto interpolate2D(FminorT const &fminor, KT const &k, int igpt, JetaT const &jeta, int jtemp,
                   int ngpt, int neta, int ntemp) {
  return fminor(0,0) * k(igpt, jeta(0)  , jtemp  ) +
         fminor(1,0) * k(igpt, jeta(0)+1, jtemp  ) +
         fminor(0,1) * k(igpt, jeta(1)  , jtemp+1) +
         fminor(1,1) * k(igpt, jeta(1)+1, jtemp+1);
}

// ----------------------------------------------------------
//
// One dimensional interpolation -- return all values along second table dimension
//
template <typename RealT, typename TableT, typename ResT>
KOKKOS_INLINE_FUNCTION
void interpolate1D(RealT val, RealT offset, RealT delta, TableT const &table,
                   ResT const &res, int tab_d1, int tab_d2) {
  RealT val0 = (val - offset) / delta;
  RealT frac = val0 - int(val0); // get fractional part
  int index = Kokkos::fmin(tab_d1-1, Kokkos::fmax(1, (int)(val0)+1)) - 1; // limit the index range
  for (int i=0; i<tab_d2; i++) {
    res(i) = table(index,i) + frac * (table(index+1,i) - table(index,i));
  }
}

template <typename FlavorT, typename PressLogT, typename TempT, typename RealT, typename VmrT, typename PlayT,
          typename TlayT, typename ColGasT, typename JtempT, typename FmajorT, typename FminorT, typename ColMixT,
          typename TropoT, typename JetaT, typename JpressT>
void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp, FlavorT const &flavor,
                   PressLogT const &press_ref_log, TempT const &temp_ref, RealT press_ref_log_delta, RealT temp_ref_min,
                   RealT temp_ref_delta, RealT press_ref_trop_log, VmrT const &vmr_ref, PlayT const &play,
                   TlayT const &tlay, ColGasT const &col_gas, JtempT const &jtemp, FmajorT const &fmajor, FminorT const &fminor,
                   ColMixT const &col_mix, TropoT const &tropo, JetaT const &jeta, JpressT const &jpress) {
  using conv::merge;

  using LayoutT = typename JtempT::array_layout;
  using DeviceT = typename JtempT::device_type;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  using pool = conv::MemPoolSingleton<RealT, LayoutT, DeviceT>;

  auto ftemp = pool::template alloc<RealT>(ncol,nlay);
  auto fpress= pool::template alloc<RealT>(ncol,nlay);

  RealT tiny = std::numeric_limits<RealT>::min();

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay, icol, ilay,
    // index and factor for temperature interpolation
    jtemp(icol,ilay) = (int) ((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta);
    jtemp(icol,ilay) = Kokkos::fmin(ntemp - 1, Kokkos::fmax(1, jtemp(icol,ilay))) - 1; // limit the index range
    ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta;

    // index and factor for pressure interpolation
    RealT locpress = 1. + (log(play(icol,ilay)) - press_ref_log(0)) / press_ref_log_delta;
    jpress(icol,ilay) = Kokkos::fmin(npres-1, Kokkos::fmax(1, (int)(locpress)));
    fpress(icol,ilay) = locpress - (RealT)(jpress(icol,ilay));

    // determine if in lower or upper part of atmosphere
    tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log;
  ));

  TIMED_KERNEL(FLATTEN_MD_KERNEL4(2,nflav,ncol,nlay, itemp, iflav, icol , ilay,
    // itropo = 0 lower atmosphere; itropo = 1 upper atmosphere
    int itropo = merge(0,1,tropo(icol,ilay));
    auto igases1 = flavor(0,iflav);
    auto igases2 = flavor(1,iflav);

    // compute interpolation fractions needed for lower, then upper reference temperature level
    // compute binary species parameter (eta) for flavor and temperature and
    //  associated interpolation index and factors
    RealT ratio_eta_half = vmr_ref(itropo,igases1+1,(jtemp(icol,ilay)+itemp)) /
                          vmr_ref(itropo,igases2+1,(jtemp(icol,ilay)+itemp));
    col_mix(itemp,iflav,icol,ilay) = col_gas(icol,ilay,igases1+1) + ratio_eta_half * col_gas(icol,ilay,igases2+1);
    RealT eta = merge(col_gas(icol,ilay,igases1+1) / col_mix(itemp,iflav,icol,ilay), 0.5,
                     col_mix(itemp,iflav,icol,ilay) > 2. * tiny);
    RealT loceta = eta * (neta-1.0);
    jeta(itemp,iflav,icol,ilay) = Kokkos::fmin((int)(loceta)+1, neta-1) - 1;
    RealT feta = fmod(loceta, 1.0);
    // compute interpolation fractions needed for minor species
    RealT ftemp_term = ((2.0 - (itemp+1)) + (2.0 * (itemp+1) - 3.0 ) * ftemp(icol,ilay));
    fminor(0,itemp,iflav,icol,ilay) = (1. - feta) * ftemp_term;
    fminor(1,itemp,iflav,icol,ilay) =       feta  * ftemp_term;
    // compute interpolation fractions needed for major species
    fmajor(0,0,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(0,itemp,iflav,icol,ilay);
    fmajor(1,0,itemp,iflav,icol,ilay) = (1. - fpress(icol,ilay)) * fminor(1,itemp,iflav,icol,ilay);
    fmajor(0,1,itemp,iflav,icol,ilay) =       fpress(icol,ilay)  * fminor(0,itemp,iflav,icol,ilay);
    fmajor(1,1,itemp,iflav,icol,ilay) =       fpress(icol,ilay)  * fminor(1,itemp,iflav,icol,ilay);
  ));

  pool::dealloc(ftemp);
  pool::dealloc(fpress);
}

template <typename TauAbsT, typename TauRayT, typename TauT, typename SsaT, typename GT>
void combine_and_reorder_2str(int ncol, int nlay, int ngpt, TauAbsT const &tau_abs, TauRayT const &tau_rayleigh,
                              TauT const &tau, SsaT const &ssa, GT const &g) {
  using RealT   = typename TauT::non_const_value_type;
  using DeviceT = typename TauT::device_type;
  using LayoutT = typename TauT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  RealT tiny = std::numeric_limits<RealT>::min();

  // int constexpr TILE_SIZE=8;
  // int colTiles = ncol / TILE_SIZE + 1;
  // int gptTiles = ngpt / TILE_SIZE + 1;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int tcol=1; tcol<=colTiles; tcol++) {
  //     for (int tgpt=1; tgpt<=gptTiles; tgpt++) {
  //       for (int itcol=1; itcol<=TILE_SIZE; itcol++) {
  //         for (int itgpt=1; itgpt<=TILE_SIZE; itgpt++) {
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<3>({ncol,nlay,ngpt}) , KOKKOS_LAMBDA (int icol, int ilay, int igpt) {
    RealT t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol);
    tau(icol,ilay,igpt) = t;
    g  (icol,ilay,igpt) = 0.;
    if(t > 2. * tiny) {
      ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t;
    } else {
      ssa(icol,ilay,igpt) = 0.;
    }
  }));
}

template <typename TlayT, typename TlevT, typename TsfcT, typename FmajorT, typename JetaT, typename TropoT,
          typename JtempT, typename JpressT, typename GpointBandsT, typename BandLimsT, typename PfracinT, typename RealT,
          typename TotplnkT, typename GpointFlavorT, typename SfcSrcT, typename LaySrcT, typename LevSrcIncT,
          typename LevSrcDecT>
void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres, int ntemp, int nPlanckTemp,
                           TlayT const &tlay, TlevT const &tlev, TsfcT const &tsfc, int sfc_lay, FmajorT const &fmajor,
                           JetaT const &jeta, TropoT const &tropo, JtempT const &jtemp, JpressT const &jpress,
                           GpointBandsT const &gpoint_bands, BandLimsT const &band_lims_gpt, PfracinT const &pfracin, RealT temp_ref_min,
                           RealT totplnk_delta, TotplnkT const &totplnk, GpointFlavorT const &gpoint_flavor, SfcSrcT const &sfc_src,
                           LaySrcT const &lay_src, LevSrcIncT const &lev_src_inc, LevSrcDecT const &lev_src_dec) {
  using conv::merge;

  using LayoutT = typename JtempT::array_layout;
  using DeviceT = typename JtempT::device_type;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;
  using pool = conv::MemPoolSingleton<RealT, LayoutT, DeviceT>;

  auto pfrac           = pool::template alloc<RealT>(ngpt,nlay,ncol);
  auto planck_function = pool::template alloc<RealT>(nbnd,nlay+1,ncol);

  // Kokkos::Array<int, 3> dims3_nlay_ncol_ngpt = {nlay,ncol,ngpt};
  // const int dims3_tot = nlay*ncol*ngpt;
  // Calculation of fraction of band's Planck irradiance associated with each g-point
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ngpt,ncol,nlay, igpt, icol, ilay,
    // itropo = 0 lower atmosphere; itropo = 1 upper atmosphere
    int itropo = merge(0,1,tropo(icol,ilay));  //WS moved itropo inside loop for GPU
    int iflav = gpoint_flavor(itropo, igpt); //eta interpolation depends on band's flavor
    // interpolation in temperature, pressure, and eta
    int jpress_loc = jpress(icol,ilay)+itropo;
    int jtemp_loc  = jtemp(icol,ilay);

    // inlining interpolate3D
    pfrac(igpt,ilay,icol) = (
       fmajor(0,0,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc  ) +
       fmajor(1,0,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc  ) +
       fmajor(0,1,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc  , jtemp_loc  ) +
       fmajor(1,1,0,iflav,icol,ilay) * pfracin(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc  ) ) +
     ( fmajor(0,0,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc+1) +
       fmajor(1,0,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc+1) +
       fmajor(0,1,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc  , jtemp_loc+1) +
       fmajor(1,1,1,iflav,icol,ilay) * pfracin(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc+1) );
  ));

  //
  // Planck function by band for the surface
  // Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
  //
  // for (int icol=1; icol<=ncol; icol++) {
  TIMED_KERNEL(Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
    interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function, Kokkos::ALL, 0 ,icol),nPlanckTemp,nbnd);
  }));
  //
  // Map to g-points
  //
  // for (int igpt=1; igpt<=ngpt; igpt++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ngpt,ncol, igpt, icol,
    sfc_src(igpt,icol) = pfrac(igpt,sfc_lay,icol) * planck_function(gpoint_bands(igpt), 0, icol);
  ));

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol,nlay, icol, ilay,
    // Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
    interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function, Kokkos::ALL,ilay,icol),nPlanckTemp,nbnd);
  ));

  //
  // Map to g-points
  //
  // Explicitly unroll a time-consuming loop here to increase instruction-level parallelism on a GPU
  // Helps to achieve higher bandwidth
  //
  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ngpt,nlay,ncol, igpt, ilay, icol,
    lay_src(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,icol);
  ));

  // compute level source irradiances for each g-point, one each for upward and downward paths
  // for (int icol=1; icol<=ncol; icol++) {
  TIMED_KERNEL(Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
    interpolate1D(tlev(icol,0), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function, Kokkos::ALL,0 ,icol),nPlanckTemp,nbnd);
  }));

  // for (int icol=1; icol<=ncol; icol++) {
  //   for (int ilay=2; ilay<=nlay+1; ilay++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol,nlay, icol, ilay,
    interpolate1D(tlev(icol,ilay+1), temp_ref_min, totplnk_delta, totplnk, Kokkos::subview(planck_function,Kokkos::ALL,ilay+1,icol),nPlanckTemp,nbnd);
  ));

  //
  // Map to g-points
  //
  // Same unrolling as mentioned before
  //
  // for (int icol=1; icol<=ncol; icol+=2) {
  //   for (int ilay=1; ilay<=nlay; ilay++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ngpt,nlay,ncol, igpt, ilay, icol,
    lev_src_dec(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay,  icol  );
    lev_src_inc(igpt,ilay,icol  ) = pfrac(igpt,ilay,icol  ) * planck_function(gpoint_bands(igpt),ilay+1,icol  );
    if (icol < ncol-1) {
      lev_src_dec(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay,  icol+1);
      lev_src_inc(igpt,ilay,icol+1) = pfrac(igpt,ilay,icol+1) * planck_function(gpoint_bands(igpt),ilay+1,icol+1);
    }
  ));

  pool::dealloc(pfrac);
  pool::dealloc(planck_function);
}

// compute Rayleigh scattering optical depths
template <typename GpointFlavorT, typename BandLimsT, typename KraylT, typename ColDryT, typename ColGasT,
          typename FminorT, typename JetaT, typename TropoT, typename JtempT, typename TauRayT>
void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp,
                          GpointFlavorT const &gpoint_flavor, BandLimsT const &band_lims_gpt, KraylT const &krayl, int idx_h2o,
                          ColDryT const &col_dry, ColGasT const &col_gas, FminorT const &fminor, JetaT const &jeta,
                          TropoT const &tropo, JtempT const &jtemp, TauRayT const &tau_rayleigh) {
  using conv::merge;
  using RealT   = typename ColDryT::non_const_value_type;
  using DeviceT = typename ColDryT::device_type;
  using LayoutT = typename ColDryT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ngpt,ncol,nlay, igpt, icol, ilay,
    int itropo = merge(0,1,tropo(icol,ilay)); // itropo = 0 lower atmosphere; itropo = 1 upper atmosphere
    int iflav = gpoint_flavor(itropo, igpt);
    // Inlining interpolate2D
    RealT k = fminor(0,0,iflav,icol,ilay) * krayl(igpt, jeta(0,iflav,icol,ilay)  , jtemp(icol,ilay)  ,itropo) +
             fminor(1,0,iflav,icol,ilay) * krayl(igpt, jeta(0,iflav,icol,ilay)+1, jtemp(icol,ilay)  ,itropo) +
             fminor(0,1,iflav,icol,ilay) * krayl(igpt, jeta(1,iflav,icol,ilay)  , jtemp(icol,ilay)+1,itropo) +
             fminor(1,1,iflav,icol,ilay) * krayl(igpt, jeta(1,iflav,icol,ilay)+1, jtemp(icol,ilay)+1,itropo);
    tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o+1)+col_dry(icol,ilay));
  ));
}

template <typename GptFlvT, typename KminorT, typename MinorLimitsT, typename MinorScalesT, typename ScaleByT,
          typename IdxMinorT, typename IdxMinorScalingT, typename KminorStartT, typename PlayT, typename TlayT,
          typename ColGasT, typename FminorT, typename JetaT, typename LayerT, typename JtempT, typename TauT>
void gas_optical_depths_minor(int max_gpt_diff, int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                              int nminor, int nminork, int idx_h2o, int idx_tropo, GptFlvT const &gpt_flv,
                              KminorT const &kminor, MinorLimitsT const &minor_limits_gpt, MinorScalesT const &minor_scales_with_density,
                              ScaleByT const &scale_by_complement, IdxMinorT const &idx_minor, IdxMinorScalingT const &idx_minor_scaling,
                              KminorStartT const &kminor_start, PlayT const &play, TlayT const &tlay, ColGasT const &col_gas,
                              FminorT const &fminor, JetaT const &jeta, LayerT const &layer_limits, JtempT const &jtemp, TauT const &tau) {
  using RealT   = typename TauT::non_const_value_type;
  using LayoutT = typename TauT::array_layout;
  using DeviceT = typename TauT::device_type;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;
  RealT constexpr PaTohPa = 0.01;

  int extent = scale_by_complement.extent(0);

  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     for (int igpt0=0; igpt0<=max_gpt_diff; igpt0++) {
  TIMED_KERNEL(Kokkos::parallel_for( mdrp_t::template get<2>({ncol,nlay}) , KOKKOS_LAMBDA (int icol, int ilay) {
    // This check skips individual columns with no pressures in range
    if ( layer_limits(icol,0) <= -1 || ilay < layer_limits(icol,0) || ilay > layer_limits(icol,1) ) {
    } else {
      RealT myplay  = play (icol,ilay);
      RealT mytlay  = tlay (icol,ilay);
      int  myjtemp = jtemp(icol,ilay);
      RealT mycol_gas_h2o = col_gas(icol,ilay,idx_h2o+1);
      RealT mycol_gas_0   = col_gas(icol,ilay,0);

// #ifndef KOKKOS_ENABLE_CUDA
      for (int imnr=0; imnr<extent; imnr++) {
// #endif
      RealT scaling = col_gas(icol,ilay,idx_minor(imnr)+1);
      if (minor_scales_with_density(imnr)) {
        // NOTE: P needed in hPa to properly handle density scaling.
        scaling = scaling * (PaTohPa * myplay/mytlay);

        if (idx_minor_scaling(imnr) > -1) {  // there is a second gas that affects this gas's absorption
          RealT mycol_gas_imnr = col_gas(icol,ilay,idx_minor_scaling(imnr)+1);
          RealT vmr_fact = 1. / mycol_gas_0;
          RealT dry_fact = 1. / (1. + mycol_gas_h2o * vmr_fact);
          // scale by density of special gas
          if (scale_by_complement(imnr)) { // scale by densities of all gases but the special one
            scaling = scaling * (1. - mycol_gas_imnr * vmr_fact * dry_fact);
          } else {
            scaling = scaling *       mycol_gas_imnr * vmr_fact * dry_fact;
          }
        }
      } // minor_scalse_with_density(imnr)

      // What is the starting point in the stored array of minor absorption coefficients?
      int minor_start = kminor_start(imnr);
      // Which gpoint range does this minor gas affect?
      int gptS = minor_limits_gpt(0,imnr);
      int gptE = minor_limits_gpt(1,imnr);
      for (int igpt=gptS; igpt<=gptE; igpt++) {
        // Interpolation of absorption coefficient and calculation of optical depth
        int iflav = gpt_flv(idx_tropo,igpt); // eta interpolation depends on flavor
        int minor_loc = minor_start + (igpt - gptS); // add offset to starting point
        // Inlined interpolate2D
        RealT kminor_loc =
          fminor(0,0,iflav,icol,ilay) * kminor(minor_loc, jeta(0,iflav,icol,ilay)  , myjtemp  ) +
          fminor(1,0,iflav,icol,ilay) * kminor(minor_loc, jeta(0,iflav,icol,ilay)+1, myjtemp  ) +
          fminor(0,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)  , myjtemp+1) +
          fminor(1,1,iflav,icol,ilay) * kminor(minor_loc, jeta(1,iflav,icol,ilay)+1, myjtemp+1);
        RealT tau_minor = kminor_loc * scaling;
// #ifdef KOKKOS_ENABLE_CUDA
//         Kokkos::atomic_add(&tau(igpt,ilay,icol), tau_minor);
// #else
        tau(igpt,ilay,icol) += tau_minor;
//#endif
      }
// #ifndef KOKKOS_ENABLE_CUDA
      }
//#endif
    }
  }));
}

// compute minor species optical depths
template <typename GpointFlavorT, typename BandLimsT, typename KmajorT, typename ColMixT, typename FmajorT,
          typename JetaT, typename TropoT, typename JtempT, typename JpressT, typename TauT>
void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, GpointFlavorT const &gpoint_flavor, BandLimsT const &band_lims_gpt, KmajorT const &kmajor,
                              ColMixT const &col_mix, FmajorT const &fmajor, JetaT const &jeta, TropoT const &tropo,
                              JtempT const &jtemp, JpressT const &jpress, TauT const &tau) {
  using conv::merge;
  using RealT   = typename TauT::non_const_value_type;
  using DeviceT = typename ColMixT::device_type;
  using LayoutT = typename ColMixT::array_layout;
  using mdrp_t  = typename conv::MDRP<LayoutT, DeviceT>;

  // optical depth calculation for major species
  // for (int ilay=1; ilay<=nlay; ilay++) {
  //   for (int icol=1; icol<=ncol; icol++) {
  //     // optical depth calculation for major species
  //     for (int igpt=1; igpt<=ngpt; igpt++) {
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ngpt,ncol,nlay, igpt, icol, ilay,
    int itropo = merge(0,1,tropo(icol,ilay));  // WS: moved inside innermost loop

    // binary species parameter (eta) and col_mix depend on band flavor
    int iflav = gpoint_flavor(itropo, igpt);
    // interpolation in temperature, pressure, and eta
    int jpress_loc = jpress(icol,ilay)+itropo;
    int jtemp_loc  = jtemp (icol,ilay);

    // inlined interpolate3D
    RealT tau_major = col_mix(0,iflav,icol,ilay) * (
      fmajor(0,0,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc  ) +
      fmajor(1,0,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc  ) +
      fmajor(0,1,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)  , jpress_loc  , jtemp_loc  ) +
      fmajor(1,1,0,iflav,icol,ilay) * kmajor(igpt, jeta(0,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc  ) ) +
      col_mix(1,iflav,icol,ilay) * (
        fmajor(0,0,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc-1, jtemp_loc+1) +
        fmajor(1,0,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc-1, jtemp_loc+1) +
        fmajor(0,1,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)  , jpress_loc  , jtemp_loc+1) +
        fmajor(1,1,1,iflav,icol,ilay) * kmajor(igpt, jeta(1,iflav,icol,ilay)+1, jpress_loc  , jtemp_loc+1) );

    tau(igpt,ilay,icol) += tau_major;
  ));
}

// Compute minor and major species opitcal depth from pre-computed interpolation coefficients
//   (jeta,jtemp,jpress)
template <typename GpointFlavorT, typename BandLimsT, typename KmajorT, typename KminorLowerT, typename KminorUpperT,
          typename MinorLimitsLowerT, typename MinorLimitsUpperT, typename MinorScalesLowerT, typename MinorScalesUpperT,
          typename ScaleByLowerT, typename ScaleByUpperT, typename IdxMinorLowerT, typename IdxMinorUpperT,
          typename IdxMinorScalingLowerT, typename IdxMinorScalingUpperT, typename KminorStartLowerT, typename KminorStartUpperT,
          typename TropoT, typename ColMixT, typename FmajorT, typename FminorT, typename PlayT, typename TlayT, typename ColGasT,
          typename JetaT, typename JtempT, typename JpressT, typename TauT>
void compute_tau_absorption(int max_gpt_diff_lower, int max_gpt_diff_upper, int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp, int nminorlower,
                            int nminorklower, int nminorupper, int nminorkupper, int idx_h2o, GpointFlavorT const &gpoint_flavor,
                            BandLimsT const &band_lims_gpt, KmajorT const &kmajor, KminorLowerT const &kminor_lower, KminorUpperT const &kminor_upper,
                            MinorLimitsLowerT const &minor_limits_gpt_lower, MinorLimitsUpperT const &minor_limits_gpt_upper, MinorScalesLowerT const &minor_scales_with_density_lower,
                            MinorScalesUpperT const &minor_scales_with_density_upper, ScaleByLowerT const &scale_by_complement_lower,
                            ScaleByUpperT const &scale_by_complement_upper, IdxMinorLowerT const &idx_minor_lower, IdxMinorUpperT const &idx_minor_upper,
                            IdxMinorScalingLowerT const &idx_minor_scaling_lower, IdxMinorScalingUpperT const &idx_minor_scaling_upper, KminorStartLowerT const &kminor_start_lower,
                            KminorStartUpperT const &kminor_start_upper, TropoT const &tropo, ColMixT const &col_mix, FmajorT const &fmajor,
                            FminorT const &fminor, PlayT const &play, TlayT const &tlay, ColGasT const &col_gas, JetaT const &jeta,
                            JtempT const &jtemp, JpressT const &jpress, TauT const &tau, bool top_at_1) {
  using RealT   = typename TauT::non_const_value_type;
  using LayoutT = typename TauT::array_layout;
  using DeviceT = typename TauT::device_type;
  using pool = conv::MemPoolSingleton<RealT, LayoutT, DeviceT>;

  auto itropo_lower= pool::template alloc<int>(ncol,2);
  auto itropo_upper= pool::template alloc<int>(ncol,2);

  int huge  = std::numeric_limits<int>::max();
  int small = std::numeric_limits<int>::min();

  if (top_at_1) {

    TIMED_KERNEL(Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
      itropo_lower(icol,1) = nlay - 1;
      {
        int minloc = huge;
        RealT mn = (RealT) huge;
        for (int i=0; i<nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = Kokkos::fmin(mn,play(icol,i));
          }
        }
        itropo_lower(icol,0) = minloc;
      }

      itropo_upper(icol,0) = 0;
      {
        int maxloc = small;
        RealT mx = (RealT) small;
        for (int i=0; i<nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = Kokkos::fmax(mx,play(icol,i));
          }
        }
        itropo_upper(icol,1) = maxloc;
      }
    }));

  } else {  // top_at_1

    TIMED_KERNEL(Kokkos::parallel_for( ncol , KOKKOS_LAMBDA ( int icol ) {
      itropo_lower(icol,0) = 0;
      {
        int minloc = huge;
        RealT mn = (RealT) huge;
        for (int i=0; i<nlay; i++) {
          if ( tropo(icol,i) ) {
            if (play(icol,i) < mn) {
              minloc = i;
            }
            mn = Kokkos::fmin(mn,play(icol,i));
          }
        }
        itropo_lower(icol,1) = minloc;
      }

      itropo_upper(icol,1) = nlay-1;
      {
        int maxloc = small;
        RealT mx = (RealT) small;
        for (int i=0; i<nlay; i++) {
          if ( !tropo(icol,i) ) {
            if (play(icol,i) > mx) {
              maxloc = i;
            }
            mx = Kokkos::fmax(mx,play(icol,i));
          }
        }
        itropo_upper(icol,0) = maxloc;
      }

    }));

  }  // top_at_1

  // ---------------------
  // Major Species
  // ---------------------
  gas_optical_depths_major(ncol, nlay, nbnd, ngpt, nflav, neta, npres, ntemp, gpoint_flavor,
                           band_lims_gpt, kmajor, col_mix, fmajor, jeta, tropo, jtemp, jpress, tau);
  // ---------------------
  // Minor Species - lower
  // ---------------------
  int idx_tropo = 0;
  gas_optical_depths_minor(max_gpt_diff_lower, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorlower, nminorklower,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_lower, minor_limits_gpt_lower,
                           minor_scales_with_density_lower, scale_by_complement_lower, idx_minor_lower,
                           idx_minor_scaling_lower, kminor_start_lower, play,  tlay, col_gas, fminor,
                           jeta, itropo_lower, jtemp, tau);
  // ---------------------
  // Minor Species - upper
  // ---------------------
  idx_tropo = 1;
  gas_optical_depths_minor(max_gpt_diff_upper, ncol, nlay, ngpt, ngas, nflav, ntemp, neta, nminorupper, nminorkupper,
                           idx_h2o, idx_tropo, gpoint_flavor, kminor_upper, minor_limits_gpt_upper,
                           minor_scales_with_density_upper, scale_by_complement_upper, idx_minor_upper,
                           idx_minor_scaling_upper, kminor_start_upper, play, tlay, col_gas, fminor,
                           jeta, itropo_upper, jtemp, tau);

  pool::dealloc(itropo_lower);
  pool::dealloc(itropo_upper);
}

#endif
