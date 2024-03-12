
#pragma once

#include "rrtmgp_const.h"

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
// Compute interpolation coefficients
// for calculations of major optical depths, minor optical depths, Rayleigh,
// and Planck fractions
void interpolation(int ncol, int nlay, int ngas, int nflav, int neta, int npres, int ntemp, int2dk const &flavor,
                   real1dk const &press_ref_log, real1dk const &temp_ref, real press_ref_log_delta, real temp_ref_min,
                   real temp_ref_delta, real press_ref_trop_log, realOff3dk const &vmr_ref, real2dk const &play,
                   real2dk const &tlay, realOff3dk const &col_gas, int2dk const &jtemp, real6dk const &fmajor, real5dk const &fminor,
                   real4dk const &col_mix, bool2dk const &tropo, int4dk const &jeta, int2dk const &jpress);

// Combine absoprtion and Rayleigh optical depths for total tau, ssa, g
void combine_and_reorder_2str(int ncol, int nlay, int ngpt, real3dk const &tau_abs, real3dk const &tau_rayleigh,
                              real3dk const &tau, real3dk const &ssa, real3dk const &g);

void compute_Planck_source(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres, int ntemp, int nPlanckTemp,
                           real2dk const &tlay, real2dk const &tlev, real1dk const &tsfc, int sfc_lay, real6dk const &fmajor,
                           int4dk const &jeta, bool2dk const &tropo, int2dk const &jtemp, int2dk const &jpress,
                           int1dk const &gpoint_bands, int2dk const &band_lims_gpt, real4dk const &pfracin, real temp_ref_min,
                           real totplnk_delta, real2dk const &totplnk, int2dk const &gpoint_flavor, real2dk const &sfc_src,
                           real3dk const &lay_src, real3dk const &lev_src_inc, real3dk const &lev_src_dec);

// compute Rayleigh scattering optical depths
void compute_tau_rayleigh(int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp,
                          int2dk const &gpoint_flavor, int2dk const &band_lims_gpt, real4dk const &krayl, int idx_h2o,
                          real2dk const &col_dry, realOff3dk const &col_gas, real5dk const &fminor, int4dk const &jeta,
                          bool2dk const &tropo, int2dk const &jtemp, real3dk const &tau_rayleigh);

void gas_optical_depths_minor(int max_gpt_diff, int ncol, int nlay, int ngpt, int ngas, int nflav, int ntemp, int neta,
                              int nminor, int nminork, int idx_h2o, int idx_tropo, int2dk const &gpt_flv,
                              real3dk const &kminor, int2dk const &minor_limits_gpt, bool1dk const &minor_scales_with_density,
                              bool1dk const &scale_by_complement, int1dk const &idx_minor, int1dk const &idx_minor_scaling,
                              int1dk const &kminor_start, real2dk const &play, real2dk const &tlay, realOff3dk const &col_gas,
                              real5dk const &fminor, int4dk const &jeta, int2dk const &layer_limits, int2dk const &jtemp, real3dk const &tau);

void gas_optical_depths_major(int ncol, int nlay, int nbnd, int ngpt, int nflav, int neta, int npres,
                              int ntemp, int2dk const &gpoint_flavor, int2dk const &band_lims_gpt, real4dk const &kmajor,
                              real4dk const &col_mix, real6dk const &fmajor, int4dk const &jeta, bool2dk const &tropo,
                              int2dk const &jtemp, int2dk const &jpress, real3dk const &tau);

void compute_tau_absorption(int max_gpt_diff_lower, int max_gpt_diff_upper, int ncol, int nlay, int nbnd, int ngpt, int ngas, int nflav, int neta, int npres, int ntemp, int nminorlower,
                            int nminorklower, int nminorupper, int nminorkupper, int idx_h2o, int2dk const &gpoint_flavor,
                            int2dk const &band_lims_gpt, real4dk const &kmajor, real3dk const &kminor_lower, real3dk const &kminor_upper,
                            int2dk const &minor_limits_gpt_lower, int2dk const &minor_limits_gpt_upper, bool1dk const &minor_scales_with_density_lower,
                            bool1dk const &minor_scales_with_density_upper, bool1dk const &scale_by_complement_lower,
                            bool1dk const &scale_by_complement_upper, int1dk const &idx_minor_lower, int1dk const &idx_minor_upper,
                            int1dk const &idx_minor_scaling_lower, int1dk const &idx_minor_scaling_upper, int1dk const &kminor_start_lower,
                            int1dk const &kminor_start_upper, bool2dk const &tropo, real4dk const &col_mix, real6dk const &fmajor,
                            real5dk const &fminor, real2dk const &play, real2dk const &tlay, realOff3dk const &col_gas, int4dk const &jeta,
                            int2dk const &jtemp, int2dk const &jpress, real3dk const &tau, bool top_at_1);

//   This function returns a single value from a subset (in gpoint) of the k table
KOKKOS_INLINE_FUNCTION
real interpolate2D(real2dk const &fminor, real3dk const &k, int igpt, int1dk const &jeta, int jtemp,
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
KOKKOS_INLINE_FUNCTION
void interpolate1D(real val, real offset, real delta, real2dk const &table,
                   real1dk const &res, int tab_d1, int tab_d2) {
  real val0 = (val - offset) / delta;
  real frac = val0 - int(val0); // get fractional part
  int index = Kokkos::fmin(tab_d1-1, Kokkos::fmax(1, (int)(val0)+1)) - 1; // limit the index range
  for (int i=0; i<tab_d2; i++) {
    res(i) = table(index,i) + frac * (table(index+1,i) - table(index,i));
  }
}


#endif
