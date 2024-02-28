
#include "mo_load_coefficients.h"
#include "YAKL_netcdf.h"

// This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
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

void load_and_init(GasOpticsRRTMGP &kdist, std::string filename, GasConcs const &available_gases) {
  yakl::SimpleNetCDF io;
  io.open(filename , yakl::NETCDF_MODE_READ);

  // Read the many arrays
  string1d     gas_names;
  string1d     gas_minor;
  string1d     identifier_minor;
  string1d     minor_gases_lower;
  string1d     minor_gases_upper;
  string1d     scaling_gas_lower;
  string1d     scaling_gas_upper;
  intHost3d    key_species;
  realHost2d   band_lims;
  intHost2d    band2gpt;
  real         press_ref_trop;
  real         temp_ref_p;
  real         temp_ref_t;
  realHost1d   press_ref;
  realHost1d   temp_ref;
  realHost3d   vmr_ref;
  realHost4d   kmajor;
  intHost2d    minor_limits_gpt_lower;
  intHost2d    minor_limits_gpt_upper;
  boolHost1d   minor_scales_with_density_lower;
  boolHost1d   minor_scales_with_density_upper;
  boolHost1d   scale_by_complement_lower;
  boolHost1d   scale_by_complement_upper;
  intHost1d    kminor_start_lower;
  intHost1d    kminor_start_upper;
  realHost3d   kminor_lower;
  realHost3d   kminor_upper;
  realHost3d   rayl_lower;
  realHost3d   rayl_upper;

  // Read in strings
  charHost2d tmp;
  tmp = charHost2d();  io.read( tmp , "gas_names"         );  gas_names         = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "gas_minor"         );  gas_minor         = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "identifier_minor"  );  identifier_minor  = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "minor_gases_lower" );  minor_gases_lower = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "minor_gases_upper" );  minor_gases_upper = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "scaling_gas_lower" );  scaling_gas_lower = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "scaling_gas_upper" );  scaling_gas_upper = char2d_to_string1d(tmp);

  io.read( key_species                     , "key_species" );
  io.read( band_lims                       , "bnd_limits_wavenumber" );
  io.read( band2gpt                        , "bnd_limits_gpt" );
  io.read( press_ref                       , "press_ref" );
  io.read( temp_ref                        , "temp_ref" );
  io.read( temp_ref_p                      , "absorption_coefficient_ref_P" );
  io.read( temp_ref_t                      , "absorption_coefficient_ref_T" );
  io.read( press_ref_trop                  , "press_ref_trop" );
  io.read( kminor_lower                    , "kminor_lower" );
  io.read( kminor_upper                    , "kminor_upper" );
  io.read( minor_limits_gpt_lower          , "minor_limits_gpt_lower" );
  io.read( minor_limits_gpt_upper          , "minor_limits_gpt_upper" );
  io.read( minor_scales_with_density_lower , "minor_scales_with_density_lower" );
  io.read( minor_scales_with_density_upper , "minor_scales_with_density_upper" );
  io.read( scale_by_complement_lower       , "scale_by_complement_lower" );
  io.read( scale_by_complement_upper       , "scale_by_complement_upper" );
  io.read( kminor_start_lower              , "kminor_start_lower" );
  io.read( kminor_start_upper              , "kminor_start_upper" );
  io.read( vmr_ref                         , "vmr_ref" );
  io.read( kmajor                          , "kmajor" );

  if (io.varExists("rayl_lower")) {
    io.read( rayl_lower , "rayl_lower" );
    io.read( rayl_upper , "rayl_upper" );
  }

  // Initialize the gas optics class with data. The calls look slightly different depending
  //   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
  // gas_optics%load() returns a string; a non-empty string indicates an error.
  if (io.varExists("totplnk")) {
    // If there's a totplnk variable in the file, then it's a longwave (internal sources) type
    realHost2d totplnk;
    realHost4d planck_frac;
    io.read( totplnk     , "totplnk"        );
    io.read( planck_frac , "plank_fraction" );
    kdist.load(available_gases, gas_names, key_species, band2gpt, band_lims, press_ref, press_ref_trop,
               temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
               gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper,
               minor_limits_gpt_lower, minor_limits_gpt_upper, minor_scales_with_density_lower,
               minor_scales_with_density_upper, scaling_gas_lower, scaling_gas_upper,
               scale_by_complement_lower, scale_by_complement_upper, kminor_start_lower,
               kminor_start_upper, totplnk, planck_frac, rayl_lower, rayl_upper);
  } else {
    // Otherwise, it's a shortwave type
    realHost1d solar_src;
    io.read( solar_src , "solar_source" );
    kdist.load(available_gases, gas_names, key_species, band2gpt, band_lims, press_ref, press_ref_trop,
               temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
               gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper,
               minor_limits_gpt_lower, minor_limits_gpt_upper, minor_scales_with_density_lower,
               minor_scales_with_density_upper, scaling_gas_lower, scaling_gas_upper,
               scale_by_complement_lower, scale_by_complement_upper, kminor_start_lower,
               kminor_start_upper, solar_src, rayl_lower, rayl_upper);
  }
  io.close();
}

#ifdef RRTMGP_ENABLE_KOKKOS
void load_and_init(GasOpticsRRTMGPK &kdist, std::string filename, GasConcsK const &available_gases) {
  conv::SimpleNetCDF io;
  io.open(filename , conv::NETCDF_MODE_READ);

  // Read the many arrays
  string1dv     gas_names;
  string1dv     gas_minor;
  string1dv     identifier_minor;
  string1dv     minor_gases_lower;
  string1dv     minor_gases_upper;
  string1dv     scaling_gas_lower;
  string1dv     scaling_gas_upper;
  intHost3dk    key_species;
  realHost2dk   band_lims;
  intHost2dk    band2gpt;
  real          press_ref_trop;
  real          temp_ref_p;
  real          temp_ref_t;
  realHost1dk   press_ref;
  realHost1dk   temp_ref;
  realHost3dk   vmr_ref;
  realHost4dk   kmajor;
  intHost2dk    minor_limits_gpt_lower;
  intHost2dk    minor_limits_gpt_upper;
  boolHost1dk   minor_scales_with_density_lower;
  boolHost1dk   minor_scales_with_density_upper;
  boolHost1dk   scale_by_complement_lower;
  boolHost1dk   scale_by_complement_upper;
  intHost1dk    kminor_start_lower;
  intHost1dk    kminor_start_upper;
  realHost3dk   kminor_lower;
  realHost3dk   kminor_upper;
  realHost3dk   rayl_lower;
  realHost3dk   rayl_upper;

  // Read in strings
  charHost2dk tmp;
  tmp = charHost2dk();  io.read( tmp , "gas_names"         );  gas_names         = char2d_to_string1d(tmp);
  tmp = charHost2dk();  io.read( tmp , "gas_minor"         );  gas_minor         = char2d_to_string1d(tmp);
  tmp = charHost2dk();  io.read( tmp , "identifier_minor"  );  identifier_minor  = char2d_to_string1d(tmp);
  tmp = charHost2dk();  io.read( tmp , "minor_gases_lower" );  minor_gases_lower = char2d_to_string1d(tmp);
  tmp = charHost2dk();  io.read( tmp , "minor_gases_upper" );  minor_gases_upper = char2d_to_string1d(tmp);
  tmp = charHost2dk();  io.read( tmp , "scaling_gas_lower" );  scaling_gas_lower = char2d_to_string1d(tmp);
  tmp = charHost2dk();  io.read( tmp , "scaling_gas_upper" );  scaling_gas_upper = char2d_to_string1d(tmp);

  io.read( key_species                     , "key_species" );
  io.read( band_lims                       , "bnd_limits_wavenumber" );
  io.read( band2gpt                        , "bnd_limits_gpt" );
  io.read( press_ref                       , "press_ref" );
  io.read( temp_ref                        , "temp_ref" );
  io.read( temp_ref_p                      , "absorption_coefficient_ref_P" );
  io.read( temp_ref_t                      , "absorption_coefficient_ref_T" );
  io.read( press_ref_trop                  , "press_ref_trop" );
  io.read( kminor_lower                    , "kminor_lower" );
  io.read( kminor_upper                    , "kminor_upper" );
  io.read( minor_limits_gpt_lower          , "minor_limits_gpt_lower" );
  io.read( minor_limits_gpt_upper          , "minor_limits_gpt_upper" );
  io.read( minor_scales_with_density_lower , "minor_scales_with_density_lower" );
  io.read( minor_scales_with_density_upper , "minor_scales_with_density_upper" );
  io.read( scale_by_complement_lower       , "scale_by_complement_lower" );
  io.read( scale_by_complement_upper       , "scale_by_complement_upper" );
  io.read( kminor_start_lower              , "kminor_start_lower" );
  io.read( kminor_start_upper              , "kminor_start_upper" );
  io.read( vmr_ref                         , "vmr_ref" );
  io.read( kmajor                          , "kmajor" );

  if (io.varExists("rayl_lower")) {
    io.read( rayl_lower , "rayl_lower" );
    io.read( rayl_upper , "rayl_upper" );
  }

  // Initialize the gas optics class with data. The calls look slightly different depending
  //   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
  // gas_optics%load() returns a string; a non-empty string indicates an error.
  if (io.varExists("totplnk")) {
    // If there's a totplnk variable in the file, then it's a longwave (internal sources) type
    realHost2dk totplnk;
    realHost4dk planck_frac;
    io.read( totplnk     , "totplnk"        );
    io.read( planck_frac , "plank_fraction" );
    kdist.load(available_gases, gas_names, key_species, band2gpt, band_lims, press_ref, press_ref_trop,
               temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
               gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper,
               minor_limits_gpt_lower, minor_limits_gpt_upper, minor_scales_with_density_lower,
               minor_scales_with_density_upper, scaling_gas_lower, scaling_gas_upper,
               scale_by_complement_lower, scale_by_complement_upper, kminor_start_lower,
               kminor_start_upper, totplnk, planck_frac, rayl_lower, rayl_upper);
  } else {
    // Otherwise, it's a shortwave type
    realHost1dk solar_src;
    io.read( solar_src , "solar_source" );
    kdist.load(available_gases, gas_names, key_species, band2gpt, band_lims, press_ref, press_ref_trop,
               temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
               gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper,
               minor_limits_gpt_lower, minor_limits_gpt_upper, minor_scales_with_density_lower,
               minor_scales_with_density_upper, scaling_gas_lower, scaling_gas_upper,
               scale_by_complement_lower, scale_by_complement_upper, kminor_start_lower,
               kminor_start_upper, solar_src, rayl_lower, rayl_upper);
  }
  io.close();
}
#endif
