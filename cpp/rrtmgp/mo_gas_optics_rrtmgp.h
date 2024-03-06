
#pragma once

#include "mo_optical_props.h"
#include "mo_source_functions.h"
#include "mo_rrtmgp_util_string.h"
#include "mo_gas_optics_kernels.h"
#include "mo_rrtmgp_constants.h"
#include "mo_rrtmgp_util_reorder.h"
#include "mo_gas_concentrations.h"

#include <iomanip>

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
//
// Class for computing spectrally-resolved gas optical properties and source functions
//   given atmopsheric physical properties (profiles of temperature, pressure, and gas concentrations)
//   The class must be initialized with data (provided as a netCDF file) before being used.
//
// Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
//   external stellar radiation (shortwave radiation in the Earth's atmosphere).
//   The variant is chosen based on what information is supplied during initialization.
//   (It might make more sense to define two sub-classes)
//
// -------------------------------------------------------------------------------------------------

#ifdef RRTMGP_ENABLE_YAKL
class GasOpticsRRTMGP : public OpticalProps {
public:
  // RRTMGP computes absorption in each band arising from
  //   two major species in each band, which are combined to make
  //     a relative mixing ratio eta and a total column amount (col_mix)
  //   contributions from zero or more minor species whose concentrations
  //     may be scaled by other components of the atmosphere
  // Absorption coefficients are interpolated from tables on a pressure/temperature/(eta) grid
  // Interpolation variables: Temperature and pressure grids
  real1d press_ref;
  real1d press_ref_log;
  real1d temp_ref;

  // Derived and stored for convenience:
  //   Min and max for temperature and pressure intepolation grids
  //   difference in ln pressure between consecutive reference levels
  //   log of reference pressure separating the lower and upper atmosphere
  real press_ref_min;
  real press_ref_max;
  real temp_ref_min;
  real temp_ref_max;
  real press_ref_log_delta;
  real temp_ref_delta;
  real press_ref_trop_log;

  int max_gpt_diff_lower;
  int max_gpt_diff_upper;

  // Major absorbers ("key species")
  //   Each unique set of major species is called a flavor.
  // Names  and reference volume mixing ratios of major gases
  string1d gas_names;  // gas names
  real3d vmr_ref;      // vmr_ref(lower or upper atmosphere, gas, temp)

  // Which two gases are in each flavor? By index
  int2d flavor;        // major species pair; (2,nflav)

  // Which flavor for each g-point? One each for lower, upper atmosphere
  int2d gpoint_flavor; // flavor = gpoint_flavor(2, g-point)

  // Major gas absorption coefficients
  real4d kmajor;       //  kmajor(g-point,eta,pressure,temperature)

  // Minor species, independently for upper and lower atmospheres
  //   Array extents in the n_minor dimension will differ between upper and lower atmospheres
  //   Each contribution has starting and ending g-points
  int2d minor_limits_gpt_lower;
  int2d minor_limits_gpt_upper;

  // Minor gas contributions might be scaled by other gas amounts; if so we need to know
  //   the total density and whether the contribution is scaled by the partner gas
  //   or its complement (i.e. all other gases)
  // Water vapor self- and foreign continua work like this, as do
  //   all collision-induced abosption pairs
  bool1d minor_scales_with_density_lower;
  bool1d minor_scales_with_density_upper;
  bool1d scale_by_complement_lower;
  bool1d scale_by_complement_upper;
  int1d idx_minor_lower;
  int1d idx_minor_upper;
  int1d idx_minor_scaling_lower;
  int1d idx_minor_scaling_upper;

  // Index into table of absorption coefficients
  int1d kminor_start_lower;
  int1d kminor_start_upper;

  // The absorption coefficients themselves
  real3d kminor_lower; // kminor_lower(n_minor,eta,temperature)
  real3d kminor_upper; // kminor_upper(n_minor,eta,temperature)

  // Rayleigh scattering coefficients
  real4d krayl; // krayl(g-point,eta,temperature,upper/lower atmosphere)

  // Planck function spectral mapping
  //   Allocated only when gas optics object is internal-source
  real4d planck_frac;   // stored fraction of Planck irradiance in band for given g-point
                        // planck_frac(g-point, eta, pressure, temperature)
  real2d totplnk;       // integrated Planck irradiance by band; (Planck temperatures,band)
  real   totplnk_delta; // temperature steps in totplnk

  // Solar source function spectral mapping
  //   Allocated only when gas optics object is external-source
  real1d solar_src; // incoming solar irradiance(g-point)

  // Ancillary
  // Index into %gas_names -- is this a key species in any band?
  bool1d is_key;

  void finalize () {
      // Free memory
      press_ref.deallocate();
      press_ref_log.deallocate();
      temp_ref.deallocate();
      gas_names.deallocate();  // gas names
      vmr_ref.deallocate();      // vmr_ref(lower or upper atmosphere, gas, temp)
      flavor.deallocate();        // major species pair; (2,nflav)
      gpoint_flavor.deallocate(); // flavor = gpoint_flavor(2, g-point)
      kmajor.deallocate();       //  kmajor(g-point,eta,pressure,temperature)
      minor_limits_gpt_lower.deallocate();
      minor_limits_gpt_upper.deallocate();
      minor_scales_with_density_lower.deallocate();
      minor_scales_with_density_upper.deallocate();
      scale_by_complement_lower.deallocate();
      scale_by_complement_upper.deallocate();
      idx_minor_lower.deallocate();
      idx_minor_upper.deallocate();
      idx_minor_scaling_lower.deallocate();
      idx_minor_scaling_upper.deallocate();
      kminor_start_lower.deallocate();
      kminor_start_upper.deallocate();
      kminor_lower.deallocate(); // kminor_lower(n_minor,eta,temperature)
      kminor_upper.deallocate(); // kminor_upper(n_minor,eta,temperature)
      krayl.deallocate(); // krayl(g-point,eta,temperature,upper/lower atmosphere)
      planck_frac.deallocate();   // stored fraction of Planck irradiance in band for given g-point
      totplnk.deallocate();       // integrated Planck irradiance by band; (Planck temperatures,band)
      solar_src.deallocate(); // incoming solar irradiance(g-point)
      is_key.deallocate();

      // Free memory in base class
      OpticalProps::finalize();
  }

  // Everything except GasConcs is on the host by default, and the available_gases.gas_name is on the host as well
  // Things will be copied to the GPU outside of this routine and stored into class variables
  void reduce_minor_arrays(GasConcs   const &available_gases,
                           string1d   const &gas_names,
                           string1d   const &gas_minor,
                           string1d   const &identifier_minor,
                           realHost3d const &kminor_atm,
                           string1d   const &minor_gases_atm,
                           intHost2d  const &minor_limits_gpt_atm,
                           boolHost1d const &minor_scales_with_density_atm,
                           string1d   const &scaling_gas_atm,
                           boolHost1d const &scale_by_complement_atm,
                           intHost1d  const &kminor_start_atm,
                           realHost3d       &kminor_atm_red,
                           string1d         &minor_gases_atm_red,
                           intHost2d        &minor_limits_gpt_atm_red,
                           boolHost1d       &minor_scales_with_density_atm_red,
                           string1d         &scaling_gas_atm_red,
                           boolHost1d       &scale_by_complement_atm_red,
                           intHost1d        &kminor_start_atm_red) {
    using yakl::intrinsics::size;

    int nm = size(minor_gases_atm,1);  // Size of the larger list of minor gases
    int tot_g = 0;
    int red_nm = 0;                    // Reduced number of minor gasses (only the ones we need)
    boolHost1d gas_is_present("gas_is_present",nm);   // Determines whether a gas in the list is needed
    // Determine the gasses needed
    for (int i=1; i <= nm; i++) {
      int idx_mnr = string_loc_in_array(minor_gases_atm(i), identifier_minor);
      gas_is_present(i) = string_in_array(gas_minor(idx_mnr),available_gases.gas_name);
      if (gas_is_present(i)) {
        tot_g = tot_g + (minor_limits_gpt_atm(2,i)-minor_limits_gpt_atm(1,i)+1);
        red_nm = red_nm + 1;
      }
    }

    // Allocate reduced arrays
    minor_gases_atm_red               = string1d  ("minor_gases_atm_red              "  ,red_nm);
    minor_scales_with_density_atm_red = boolHost1d("minor_scales_with_density_atm_red"  ,red_nm);
    scaling_gas_atm_red               = string1d  ("scaling_gas_atm_red              "  ,red_nm);
    scale_by_complement_atm_red       = boolHost1d("scale_by_complement_atm_red      "  ,red_nm);
    kminor_start_atm_red              = intHost1d ("kminor_start_atm_red             "  ,red_nm);
    minor_limits_gpt_atm_red          = intHost2d ("minor_limits_gpt_atm_red         ",2,red_nm);
    kminor_atm_red                    = realHost3d("kminor_atm_red                   ",tot_g , size(kminor_atm,2) , size(kminor_atm,3));

    if (red_nm == nm) {
      // If the gasses listed exactly matches the gasses needed, just copy it
      minor_gases_atm              .deep_copy_to(minor_gases_atm_red              );
      scaling_gas_atm              .deep_copy_to(scaling_gas_atm_red              );
      kminor_atm                   .deep_copy_to(kminor_atm_red                   );
      minor_limits_gpt_atm         .deep_copy_to(minor_limits_gpt_atm_red         );
      minor_scales_with_density_atm.deep_copy_to(minor_scales_with_density_atm_red);
      scale_by_complement_atm      .deep_copy_to(scale_by_complement_atm_red      );
      kminor_start_atm             .deep_copy_to(kminor_start_atm_red             );
    } else {
      // Otherwise, pack into reduced arrays
      int slot = 1;
      for (int i=1; i <= nm; i++) {
        if (gas_is_present(i)) {
          minor_gases_atm_red              (slot) = minor_gases_atm              (i);
          scaling_gas_atm_red              (slot) = scaling_gas_atm              (i);
          minor_scales_with_density_atm_red(slot) = minor_scales_with_density_atm(i);
          scale_by_complement_atm_red      (slot) = scale_by_complement_atm      (i);
          kminor_start_atm_red             (slot) = kminor_start_atm             (i);
          slot++;
        }
      }

      slot = 0;
      int n_elim = 0;
      for (int i=1 ; i <= nm ; i++) {
        int ng = minor_limits_gpt_atm(2,i)-minor_limits_gpt_atm(1,i)+1;
        if (gas_is_present(i)) {
          slot = slot + 1;
          minor_limits_gpt_atm_red   (1,slot) = minor_limits_gpt_atm(1,i);
          minor_limits_gpt_atm_red   (2,slot) = minor_limits_gpt_atm(2,i);
          kminor_start_atm_red         (slot) = kminor_start_atm(i)-n_elim;
          for (int l=1 ; l <= size(kminor_atm,3) ; l++ ) {
            for (int k=1 ; k <= size(kminor_atm,2) ; k++ ) {
              for (int j=1 ; j <= ng ; j++) {
                kminor_atm_red(kminor_start_atm_red(slot)+j-1,k,l) = kminor_atm(kminor_start_atm(i)+j-1,k,l);
              }
            }
          }
        } else {
          n_elim = n_elim + ng;
        }
      }
    }
  }



  // create index list for extracting col_gas needed for minor gas optical depth calculations
  void create_idx_minor(string1d const &gas_names, string1d const &gas_minor, string1d const &identifier_minor,
                        string1d const &minor_gases_atm, intHost1d &idx_minor_atm) {
    using yakl::intrinsics::size;

    idx_minor_atm = intHost1d("idx_minor_atm",size(minor_gases_atm,1));
    for (int imnr=1 ; imnr <= size(minor_gases_atm,1) ; imnr++) {
      // Find identifying string for minor species in list of possible identifiers (e.g. h2o_slf)
      int idx_mnr     = string_loc_in_array(minor_gases_atm(imnr), identifier_minor);
      // Find name of gas associated with minor species identifier (e.g. h2o)
      idx_minor_atm(imnr) = string_loc_in_array(gas_minor(idx_mnr), gas_names);
    }
  }



  // create index for special treatment in density scaling of minor gases
  void create_idx_minor_scaling(string1d const &gas_names, string1d const &scaling_gas_atm,
                                intHost1d &idx_minor_scaling_atm) {
    using yakl::intrinsics::size;

    idx_minor_scaling_atm = intHost1d("idx_minor_scaling_atm",size(scaling_gas_atm,1));
    for (int imnr=1 ; imnr <= size(scaling_gas_atm,1) ; imnr++) {
      // This will be -1 if there's no interacting gas
      idx_minor_scaling_atm(imnr) = string_loc_in_array(scaling_gas_atm(imnr), gas_names);
    }
  }



  void create_key_species_reduce(string1d const &gas_names, string1d const &gas_names_red, intHost3d const &key_species,
                                 intHost3d &key_species_red, boolHost1d &key_species_present_init) {
    using yakl::intrinsics::size;

    int np = size(key_species,1);
    int na = size(key_species,2);
    int nt = size(key_species,3);
    key_species_red = intHost3d("key_species_red",np,na,nt);
    key_species_present_init = boolHost1d("key_species_present_init",size(gas_names,1));

    for (int i=1 ; i <= size(gas_names,1) ; i++) {
      key_species_present_init(i) = true;
    }

    for (int ip=1 ; ip <= np ; ip++) {
      for (int ia=1 ; ia <= na ; ia++) {
        for (int it=1 ; it <= nt ; it++) {
          if (key_species(ip,ia,it) != 0) {
            key_species_red(ip,ia,it) = string_loc_in_array(gas_names(key_species(ip,ia,it)),gas_names_red);
            if (key_species_red(ip,ia,it) == -1) {
              key_species_present_init(key_species(ip,ia,it)) = false;
            }
          } else {
            key_species_red(ip,ia,it) = 0;
          }
        }
      }
    }
  }



  // Create flavor list
  // An unordered array of extent (2,:) containing all possible pairs of key species used in either upper or lower atmos
  void create_flavor(intHost3d const &key_species, intHost2d &flavor) {
    using yakl::intrinsics::size;

    // prepare list of key_species
    int i = 1;
    intHost2d key_species_list("key_species_list",2,size(key_species,3)*2);
    for (int ibnd=1 ; ibnd <= size(key_species,3) ; ibnd++) {
      for (int iatm=1 ; iatm <= size(key_species,1) ; iatm++) {
        key_species_list(1,i) = key_species(1,iatm,ibnd);
        key_species_list(2,i) = key_species(2,iatm,ibnd);
        i = i + 1;
      }
    }
    // rewrite single key_species pairs
    for (int i=1 ; i <= size(key_species_list,2) ; i++) {
      if (key_species_list(1,i) == 0 && key_species_list(2,i) == 0) {
        key_species_list(1,i) = 2;
        key_species_list(2,i) = 2;
      }
    }
    // count unique key species pairs
    int iflavor = 0;
    for (int i=1; i <= size(key_species_list,2) ; i++) {
      // Loop through previous pairs. Only increment iflavor if we haven't seen this pair before
      bool unique = true;
      for (int j=1; j <= i-1 ; j++) {
        if ( key_species_list(1,j) == key_species_list(1,i) && key_species_list(2,j) == key_species_list(2,i) ) {
          unique = false;
        }
      }
      if (unique) { iflavor = iflavor + 1; }
    }
    // fill flavors
    flavor = intHost2d("flavor",2,iflavor);
    iflavor = 0;
    for (int i=1 ; i <= size(key_species_list,2) ; i++) {
      bool unique = true;
      for (int j=1; j <= i-1 ; j++) {
        if ( key_species_list(1,j) == key_species_list(1,i) && key_species_list(2,j) == key_species_list(2,i) ) {
          unique = false;
        }
      }
      if (unique) {
        iflavor = iflavor + 1;
        flavor(1,iflavor) = key_species_list(1,i);
        flavor(2,iflavor) = key_species_list(2,i);
      }
    }
  }



  // create gpoint_flavor list: a map pointing from each g-point to the corresponding entry in the "flavor list"
  void create_gpoint_flavor(intHost3d const &key_species, intHost1d const &gpt2band, intHost2d const &flavor,
                            intHost2d &gpoint_flavor) {
    using yakl::intrinsics::size;

    int ngpt = size(gpt2band,1);
    gpoint_flavor = intHost2d("gpoint_flavor",2,ngpt);
    for (int igpt=1 ; igpt <= ngpt ; igpt++) {
      for (int iatm = 1 ; iatm <= 2 ; iatm++) {
        int key_species_pair2flavor = 0;
        for (int iflav=1 ; iflav <= size(flavor,2) ; iflav++) {
          int key_species1 = key_species(1,iatm,gpt2band(igpt));
          int key_species2 = key_species(2,iatm,gpt2band(igpt));
          if (key_species1 == 0 && key_species2 == 0) {
            key_species1 = 2;
            key_species2 = 2;
          }
          if ( flavor(1,iflav) == key_species1 && flavor(2,iflav) == key_species2 ) {
            key_species_pair2flavor = iflav;
          }
        }
        gpoint_flavor(iatm,igpt) = key_species_pair2flavor;
      }
    }
  }



  // Initialize absorption coefficient arrays,
  //   including Rayleigh scattering tables if provided (allocated)
  void init_abs_coeffs(GasConcs   const &available_gases,
                       string1d   const &gas_names,
                       intHost3d  const &key_species,
                       intHost2d  const &band2gpt,
                       realHost2d const &band_lims_wavenum,
                       realHost1d const &press_ref,
                       realHost1d const &temp_ref,
                       real              press_ref_trop,
                       real              temp_ref_p,
                       real              temp_ref_t,
                       realHost3d const &vmr_ref,
                       realHost4d const &kmajor,
                       realHost3d const &kminor_lower,
                       realHost3d const &kminor_upper,
                       string1d   const &gas_minor,
                       string1d   const &identifier_minor,
                       string1d   const &minor_gases_lower,
                       string1d   const &minor_gases_upper,
                       intHost2d  const &minor_limits_gpt_lower,
                       intHost2d  const &minor_limits_gpt_upper,
                       boolHost1d const &minor_scales_with_density_lower,
                       boolHost1d const &minor_scales_with_density_upper,
                       string1d   const &scaling_gas_lower,
                       string1d   const &scaling_gas_upper,
                       boolHost1d const &scale_by_complement_lower,
                       boolHost1d const &scale_by_complement_upper,
                       intHost1d  const &kminor_start_lower,
                       intHost1d  const &kminor_start_upper,
                       realHost3d const &rayl_lower,
                       realHost3d const &rayl_upper) {
    using yakl::intrinsics::count;
    using yakl::intrinsics::pack;
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    OpticalProps::init(band_lims_wavenum.createDeviceCopy(), band2gpt.createDeviceCopy());

    // Which gases known to the gas optics are present in the host model (available_gases)?
    int ngas = size(gas_names,1);
    boolHost1d gas_is_present("gas_is_present",ngas);
    for (int i=1 ; i <= ngas ; i++) {
      gas_is_present(i) = string_in_array(gas_names(i), available_gases.gas_name);
    }
    // Now the number of gases is the union of those known to the k-distribution and provided by the host model
    ngas = count(gas_is_present);

    // Initialize the gas optics object, keeping only those gases known to the gas optics and also present in the host model
    this->gas_names = pack(gas_names,gas_is_present);

    realHost3d vmr_ref_red("vmr_ref_red",size(vmr_ref,1),{0,ngas}, size(vmr_ref,3));
    // Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    for (int k=1 ; k <= size(vmr_ref,3) ; k++) {
      for (int j=1 ; j <= size(vmr_ref,1) ; j++) {
        vmr_ref_red(j,0,k) = vmr_ref(j,1,k);
      }
    }
    for (int i=1 ; i <= ngas ; i++) {
      int idx = string_loc_in_array(this->gas_names(i), gas_names);
      for (int k=1 ; k <= size(vmr_ref,3) ; k++) {
        for (int j=1 ; j <= size(vmr_ref,1) ; j++) {
          vmr_ref_red(j,i,k) = vmr_ref(j,idx+1,k);
        }
      }
    }
    // Allocate class copy, and deep copy to the class data member
    this->vmr_ref = real3d("vmr_ref",size(vmr_ref,1),{0,ngas}, size(vmr_ref,3));
    vmr_ref_red.deep_copy_to(this->vmr_ref);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // REDUCE MINOR ARRAYS SO VARIABLES ONLY CONTAIN MINOR GASES THAT ARE AVAILABLE
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // LOWER MINOR GASSES
    string1d minor_gases_lower_red;
    string1d scaling_gas_lower_red;
    realHost3d kminor_lower_red;
    intHost2d  minor_limits_gpt_lower_red;
    boolHost1d minor_scales_with_density_lower_red;
    boolHost1d scale_by_complement_lower_red;
    intHost1d  kminor_start_lower_red;

    reduce_minor_arrays(available_gases, gas_names, gas_minor, identifier_minor, kminor_lower, minor_gases_lower,
                        minor_limits_gpt_lower, minor_scales_with_density_lower, scaling_gas_lower,
                        scale_by_complement_lower, kminor_start_lower, kminor_lower_red, minor_gases_lower_red,
                        minor_limits_gpt_lower_red, minor_scales_with_density_lower_red, scaling_gas_lower_red,
                        scale_by_complement_lower_red, kminor_start_lower_red);

    this->kminor_lower                    = kminor_lower_red                   .createDeviceCopy();
    this->minor_limits_gpt_lower          = minor_limits_gpt_lower_red         .createDeviceCopy();
    this->minor_scales_with_density_lower = minor_scales_with_density_lower_red.createDeviceCopy();
    this->scale_by_complement_lower       = scale_by_complement_lower_red      .createDeviceCopy();
    this->kminor_start_lower              = kminor_start_lower_red             .createDeviceCopy();

    // Find the largest number of g-points per band
    this->max_gpt_diff_lower = minor_limits_gpt_lower_red(2,1) - minor_limits_gpt_lower_red(1,1);
    for (int i=2; i<=size(minor_limits_gpt_lower_red,2); i++) {
      this->max_gpt_diff_lower = std::max( this->max_gpt_diff_lower , minor_limits_gpt_lower_red(2,i) - minor_limits_gpt_lower_red(1,i) );
    }

    // UPPER MINOR GASSES
    string1d minor_gases_upper_red;
    string1d scaling_gas_upper_red;
    realHost3d kminor_upper_red;
    intHost2d  minor_limits_gpt_upper_red;
    boolHost1d minor_scales_with_density_upper_red;
    boolHost1d scale_by_complement_upper_red;
    intHost1d  kminor_start_upper_red;

    reduce_minor_arrays(available_gases, gas_names, gas_minor, identifier_minor, kminor_upper, minor_gases_upper,
                        minor_limits_gpt_upper, minor_scales_with_density_upper, scaling_gas_upper,
                        scale_by_complement_upper, kminor_start_upper, kminor_upper_red, minor_gases_upper_red,
                        minor_limits_gpt_upper_red, minor_scales_with_density_upper_red, scaling_gas_upper_red,
                        scale_by_complement_upper_red, kminor_start_upper_red);

    this->kminor_upper                    = kminor_upper_red                   .createDeviceCopy();
    this->minor_limits_gpt_upper          = minor_limits_gpt_upper_red         .createDeviceCopy();
    this->minor_scales_with_density_upper = minor_scales_with_density_upper_red.createDeviceCopy();
    this->scale_by_complement_upper       = scale_by_complement_upper_red      .createDeviceCopy();
    this->kminor_start_upper              = kminor_start_upper_red             .createDeviceCopy();

    // Find the largest number of g-points per band
    this->max_gpt_diff_upper = minor_limits_gpt_upper_red(2,1) - minor_limits_gpt_upper_red(1,1);
    for (int i=2; i<=size(minor_limits_gpt_upper_red,2); i++) {
      this->max_gpt_diff_upper = std::max( this->max_gpt_diff_upper , minor_limits_gpt_upper_red(2,i) - minor_limits_gpt_upper_red(1,i) );
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // HANDLE ARRAYS NOT REDUCED BY THE PRESENCE, OR LACK THEREOF, OF A GAS
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    this->press_ref = real1d("press_ref",size(press_ref,1));
    this->temp_ref  = real1d("temp_ref" ,size(temp_ref ,1));
    this->kmajor    = real4d("kmajor"   ,size(kmajor,1),size(kmajor,2),size(kmajor,3),size(kmajor,4));
    press_ref.deep_copy_to(this->press_ref);
    temp_ref .deep_copy_to(this->temp_ref );
    kmajor   .deep_copy_to(this->kmajor   );

    // Process rayl_lower and rayl_upper into a combined this->krayl
    if (allocated(rayl_lower) != allocated(rayl_upper)) {
      stoprun("rayl_lower and rayl_upper must have the same allocation status");
    }
    if (allocated(rayl_lower)) {
      realHost4d krayltmp("krayltmp",size(rayl_lower,1),size(rayl_lower,2),size(rayl_lower,3),2);
      for (int k=1 ; k <= size(rayl_lower,3) ; k++ ) {
        for (int j=1 ; j <= size(rayl_lower,2) ; j++ ) {
          for (int i=1 ; i <= size(rayl_lower,1) ; i++ ) {
            krayltmp(i,j,k,1) = rayl_lower(i,j,k);
            krayltmp(i,j,k,2) = rayl_upper(i,j,k);
          }
        }
      }
      this->krayl = krayltmp.createDeviceCopy();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // POST PROCESSING
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // creates log reference pressure
    this->press_ref_log = real1d("press_ref_log",size(this->press_ref,1));
    YAKL_SCOPE( press_ref_loc     , this->press_ref     );
    YAKL_SCOPE( press_ref_log_loc , this->press_ref_log );
    // Running a kernel because it's more convenient in this case
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<1>( size(this->press_ref,1) ) , YAKL_LAMBDA (int i) {
      press_ref_log_loc(i) = log(press_ref_loc(i));
    });

    // log scale of reference pressure (this is a scalar, not an array)
    this->press_ref_trop_log = log(press_ref_trop);

    // Get index of gas (if present) for determining col_gas
    intHost1d idx_minor_lower_tmp;
    intHost1d idx_minor_upper_tmp;
    create_idx_minor(this->gas_names, gas_minor, identifier_minor, minor_gases_lower_red, idx_minor_lower_tmp);
    create_idx_minor(this->gas_names, gas_minor, identifier_minor, minor_gases_upper_red, idx_minor_upper_tmp);
    this->idx_minor_lower = idx_minor_lower_tmp.createDeviceCopy();
    this->idx_minor_upper = idx_minor_upper_tmp.createDeviceCopy();
    // Get index of gas (if present) that has special treatment in density scaling
    intHost1d idx_minor_scaling_lower_tmp;
    intHost1d idx_minor_scaling_upper_tmp;
    create_idx_minor_scaling(this->gas_names, scaling_gas_lower_red, idx_minor_scaling_lower_tmp);
    create_idx_minor_scaling(this->gas_names, scaling_gas_upper_red, idx_minor_scaling_upper_tmp);
    this->idx_minor_scaling_lower = idx_minor_scaling_lower_tmp.createDeviceCopy();
    this->idx_minor_scaling_upper = idx_minor_scaling_upper_tmp.createDeviceCopy();

    // create flavor list
    // Reduce (remap) key_species list; checks that all key gases are present in incoming
    boolHost1d key_species_present_init;
    intHost3d  key_species_red;
    create_key_species_reduce(gas_names, this->gas_names, key_species, key_species_red, key_species_present_init);
    // create flavor and gpoint_flavor lists
    intHost2d flavor_tmp;
    intHost2d gpoint_flavor_tmp;
    create_flavor       (key_species_red, flavor_tmp);
    create_gpoint_flavor(key_species_red, this->get_gpoint_bands().createHostCopy(), flavor_tmp, gpoint_flavor_tmp);
    this->flavor        = int2d("flavor",size(       flavor_tmp,1),size(       flavor_tmp,2));
    this->gpoint_flavor = int2d("flavor",size(gpoint_flavor_tmp,1),size(gpoint_flavor_tmp,2));
    flavor_tmp       .deep_copy_to(this->flavor       );
    gpoint_flavor_tmp.deep_copy_to(this->gpoint_flavor);

    // minimum, maximum reference temperature, pressure -- assumes low-to-high ordering
    //   for T, high-to-low ordering for p
    this->temp_ref_min  = temp_ref (1);
    this->temp_ref_max  = temp_ref (size(temp_ref ,1));
    this->press_ref_min = press_ref(size(press_ref,1));
    this->press_ref_max = press_ref(1);

    // creates press_ref_log, temp_ref_delta
    this->press_ref_log_delta = (log(this->press_ref_min)-log(this->press_ref_max))/(size(press_ref,1)-1);
    this->temp_ref_delta      = (this->temp_ref_max-this->temp_ref_min)/(size(temp_ref,1)-1);

    // Which species are key in one or more bands?
    //   this->flavor is an index into this->gas_names
    this->is_key = bool1d("is_key",this->get_ngas());
    YAKL_SCOPE( is_key_loc , this->is_key );
    YAKL_SCOPE( flavor_loc , this->flavor );
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<1>( this->get_ngas() ) , YAKL_LAMBDA (int i) {
      is_key_loc(i) = false;
    });
    // do j = 1, size(this%flavor, 2)
    //   do i = 1, size(this%flavor, 1) ! extents should be 2
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>( size(this->flavor,2) , size(this->flavor,1) ) , YAKL_LAMBDA (int j, int i) {
      if (flavor_loc(i,j) != 0) { is_key_loc(flavor_loc(i,j)) = true; }
    });
  }



  // Initialize object based on data read from netCDF file however the user desires.
  //  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  // This interface is for the internal-sources object -- includes Plank functions and fractions
  void load(GasConcs   const &available_gases,
            string1d   const &gas_names,
            intHost3d  const &key_species,
            intHost2d  const &band2gpt,
            realHost2d const &band_lims_wavenum,
            realHost1d const &press_ref,
            real              press_ref_trop,
            realHost1d const &temp_ref,
            real              temp_ref_p,
            real              temp_ref_t,
            realHost3d const &vmr_ref,
            realHost4d const &kmajor,
            realHost3d const &kminor_lower,
            realHost3d const &kminor_upper,
            string1d   const &gas_minor,
            string1d   const &identifier_minor,
            string1d   const &minor_gases_lower,
            string1d   const &minor_gases_upper,
            intHost2d  const &minor_limits_gpt_lower,
            intHost2d  const &minor_limits_gpt_upper,
            boolHost1d const &minor_scales_with_density_lower,
            boolHost1d const &minor_scales_with_density_upper,
            string1d   const &scaling_gas_lower,
            string1d   const &scaling_gas_upper,
            boolHost1d const &scale_by_complement_lower,
            boolHost1d const &scale_by_complement_upper,
            intHost1d  const &kminor_start_lower,
            intHost1d  const &kminor_start_upper,
            realHost2d const &totplnk,
            realHost4d const &planck_frac,
            realHost3d const &rayl_lower,
            realHost3d const &rayl_upper) {
    using yakl::intrinsics::size;

    init_abs_coeffs(available_gases, gas_names, key_species, band2gpt, band_lims_wavenum, press_ref, temp_ref,
                    press_ref_trop, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
                    gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,
                    minor_limits_gpt_upper, minor_scales_with_density_lower, minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,
                    kminor_start_lower, kminor_start_upper, rayl_lower, rayl_upper);

    // Planck function tables
    this->totplnk = real2d("totplnk",size(totplnk,1),size(totplnk,2));
    this->planck_frac = real4d("planck_frac",size(planck_frac,1),size(planck_frac,2),size(planck_frac,3),size(planck_frac,4));
    totplnk    .deep_copy_to(this->totplnk);
    planck_frac.deep_copy_to(this->planck_frac);
    // Temperature steps for Planck function interpolation
    //   Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
    //   Planck grid and the Planck grid is equally spaced
    this->totplnk_delta = (this->temp_ref_max - this->temp_ref_min) / (size(this->totplnk,1)-1);
  }



  // Initialize object based on data read from netCDF file however the user desires.
  //  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  // This interface is for the external-sources object -- includes TOA source function table
  void load(GasConcs   const &available_gases,
            string1d   const &gas_names,
            intHost3d  const &key_species,
            intHost2d  const &band2gpt,
            realHost2d const &band_lims_wavenum,
            realHost1d const &press_ref,
            real              press_ref_trop,
            realHost1d const &temp_ref,
            real              temp_ref_p,
            real              temp_ref_t,
            realHost3d const &vmr_ref,
            realHost4d const &kmajor,
            realHost3d const &kminor_lower,
            realHost3d const &kminor_upper,
            string1d   const &gas_minor,
            string1d   const &identifier_minor,
            string1d   const &minor_gases_lower,
            string1d   const &minor_gases_upper,
            intHost2d  const &minor_limits_gpt_lower,
            intHost2d  const &minor_limits_gpt_upper,
            boolHost1d const &minor_scales_with_density_lower,
            boolHost1d const &minor_scales_with_density_upper,
            string1d   const &scaling_gas_lower,
            string1d   const &scaling_gas_upper,
            boolHost1d const &scale_by_complement_lower,
            boolHost1d const &scale_by_complement_upper,
            intHost1d  const &kminor_start_lower,
            intHost1d  const &kminor_start_upper,
            realHost1d const &solar_src,
            realHost3d const &rayl_lower,
            realHost3d const &rayl_upper) {
    using yakl::intrinsics::size;

    init_abs_coeffs(available_gases,  gas_names, key_species, band2gpt, band_lims_wavenum, press_ref, temp_ref,
                    press_ref_trop, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
                    gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,
                    minor_limits_gpt_upper, minor_scales_with_density_lower, minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,
                    kminor_start_lower, kminor_start_upper, rayl_lower, rayl_upper);

    // Solar source table init
    this->solar_src = real1d("solar_src",size(solar_src,1));
    solar_src.deep_copy_to(this->solar_src);
    this->totplnk_delta = 0.;
  }



  // Two functions to define array sizes needed by gas_optics()
  int get_ngas() const { return yakl::intrinsics::size(this->gas_names,1); }



  // return the number of distinct major gas pairs in the spectral bands (referred to as
  // "flavors" - all bands have a flavor even if there is one or no major gas)
  int get_nflav() const { return yakl::intrinsics::size(this->flavor,2); }



  string1d get_gases() const { return this->gas_names; }



  // return the minimum pressure on the interpolation grids
  real get_press_min() const { return this->press_ref_min; }



  // return the maximum pressure on the interpolation grids
  real get_press_max() const { return this->press_ref_max; }



  // return the minimum temparature on the interpolation grids
  real get_temp_min()const  { return this->temp_ref_min; }



  // return the maximum temparature on the interpolation grids
  real get_temp_max() const { return this->temp_ref_max; }



  int get_neta() const { return yakl::intrinsics::size(this->kmajor,2); }



  // return the number of pressures in reference profile
  //   absorption coefficient table is one bigger since a pressure is repeated in upper/lower atmos
  int get_npres() const { return yakl::intrinsics::size(this->kmajor,3)-1; }



  int get_ntemp() const { return yakl::intrinsics::size(this->kmajor,4); }



  // return the number of temperatures for Planck function
  int get_nPlanckTemp() const { return yakl::intrinsics::size(this->totplnk,1); }



  // Function to define names of key and minor gases to be used by gas_optics().
  // The final list gases includes those that are defined in gas_optics_specification
  // and are provided in ty_gas_concs.
  string1d get_minor_list(GasConcs const &gas_desc, int ngas, string1d const &name_spec) const {
    using yakl::intrinsics::pack;
    using yakl::intrinsics::size;
    // List of minor gases to be used in gas_optics()
    boolHost1d gas_is_present("gas_is_present",size(name_spec,1));
    for (int igas=1 ; igas <= this->get_ngas() ; igas++) {
      gas_is_present(igas) = string_in_array(name_spec(igas), gas_desc.gas_name);
    }
    return pack(this->gas_names, gas_is_present);
  }



  // return true if initialized for internal sources, false otherwise
  bool source_is_internal() const { return yakl::intrinsics::allocated(this->totplnk) && yakl::intrinsics::allocated(this->planck_frac); }



  // return true if initialized for external sources, false otherwise
  bool source_is_external() const { return yakl::intrinsics::allocated(this->solar_src); }



  // Ensure that every key gas required by the k-distribution is present in the gas concentration object
  void check_key_species_present(GasConcs const &gas_desc) const {
    using yakl::intrinsics::pack;
    using yakl::intrinsics::size;

    string1d key_gas_names = pack(this->gas_names, this->is_key.createHostCopy());
    for (int igas=1 ; igas <= size(key_gas_names,1) ; igas++) {
      if (! string_in_array(key_gas_names(igas), gas_desc.gas_name)) {
        stoprun("gas required by k-distribution is not present in the GasConcs object");
      }
    }
  }



  // Compute gas optical depth and Planck source functions, given temperature, pressure, and composition
  template <class T>
  void gas_optics(const int ncol, const int nlay,
                  bool top_at_1, real2d const &play, real2d const &plev, real2d const &tlay, real1d const &tsfc,
                  GasConcs const &gas_desc, T &optical_props, SourceFuncLW &sources,
                  real2d const &col_dry=real2d(), real2d const &tlev=real2d()) {
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;
    using yakl::intrinsics::any;
    using yakl::componentwise::operator>;
    using yakl::componentwise::operator<;

    int ngpt  = this->get_ngpt();
    int nband = this->get_nband();
    // Interpolation coefficients for use in source function
    int2d  jtemp ("jtemp"                         ,ncol,nlay);
    int2d  jpress("jpress"                        ,ncol,nlay);
    bool2d tropo ("tropo"                         ,ncol,nlay);
    real6d fmajor("fmajor",2,2,2,this->get_nflav(),ncol,nlay);
    int4d  jeta  ("jeta"  ,2    ,this->get_nflav(),ncol,nlay);
    // Gas optics
    compute_gas_taus(top_at_1, ncol, nlay, ngpt, nband, play, plev, tlay, gas_desc, optical_props, jtemp, jpress,
                     jeta, tropo, fmajor, col_dry);

    // External source -- check arrays sizes and values
    // input data sizes and values
    if (size(tsfc,1) != ncol) { stoprun("gas_optics(): array tsfc has wrong size"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(tsfc < this->temp_ref_min) || any(tsfc > this->temp_ref_max)) {
        stoprun("gas_optics(): array tsfc has values outside range");
      }
    #endif

    if (allocated(tlev)) {
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (any(tlev < this->temp_ref_min) || any(tlev > this->temp_ref_max)) {
          stoprun("gas_optics(): array tlev has values outside range");
        }
      #endif
    }

    // output extents
    if (sources.get_ncol() != ncol || sources.get_nlay() != nlay || sources.get_ngpt() != ngpt) {
      stoprun("gas_optics%gas_optics: source function arrays inconsistently sized");
    }

    // Interpolate source function
    this->source(top_at_1, ncol, nlay, nband, ngpt, play, plev, tlay, tsfc, jtemp, jpress, jeta, tropo, fmajor, sources, tlev);
  }



  // Compute gas optical depth given temperature, pressure, and composition
  template <class T>
  void gas_optics(const int ncol, const int nlay,
                  bool top_at_1, real2d const &play, real2d const &plev, real2d const &tlay, GasConcs const &gas_desc,
                  T &optical_props, real2d &toa_src, real2d const &col_dry=real2d()) {
    using yakl::intrinsics::size;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    int ngpt  = this->get_ngpt();
    int nband = this->get_nband();
    int ngas  = this->get_ngas();
    int nflav = get_nflav();

    // Interpolation coefficients for use in source function
    int2d  jtemp ("jtemp"                         ,ncol,nlay);
    int2d  jpress("jpress"                        ,ncol,nlay);
    bool2d tropo ("tropo"                         ,ncol,nlay);
    real6d fmajor("fmajor",2,2,2,this->get_nflav(),ncol,nlay);
    int4d  jeta  ("jeta  ",2    ,this->get_nflav(),ncol,nlay);
    // Gas optics
    compute_gas_taus(top_at_1, ncol, nlay, ngpt, nband, play, plev, tlay, gas_desc, optical_props, jtemp, jpress, jeta,
                     tropo, fmajor, col_dry);

    // External source function is constant
    if (size(toa_src,1) != ncol || size(toa_src,2) != ngpt) { stoprun("gas_optics(): array toa_src has wrong size"); }

    YAKL_SCOPE( solar_src_loc , this->solar_src );
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      toa_src(icol,igpt) = solar_src_loc(igpt);
    });
  }



  // Returns optical properties and interpolation coefficients
  template <class T>
  void compute_gas_taus(bool top_at_1, int ncol, int nlay, int ngpt, int nband, real2d const &play, real2d const &plev, real2d const &tlay,
                        GasConcs const &gas_desc, T &optical_props, int2d const &jtemp, int2d const &jpress, int4d const &jeta,
                        bool2d const &tropo, real6d const &fmajor, real2d const &col_dry=real2d() ) {
    using yakl::intrinsics::lbound;
    using yakl::intrinsics::ubound;
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;
    using yakl::intrinsics::any;
    using yakl::componentwise::operator>;
    using yakl::componentwise::operator<;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;
    using yakl::COLON;

    // Number of molecules per cm^2
    real3d tau         ("tau"         ,ngpt,nlay,ncol);
    real3d tau_rayleigh("tau_rayleigh",ngpt,nlay,ncol);
    // Interpolation variables used in major gas but not elsewhere, so don't need exporting
    real3d vmr         ("vmr"         ,ncol,nlay,this->get_ngas());
    real3d col_gas     ("col_gas"     ,ncol,nlay,{0,this->get_ngas()});
    real4d col_mix     ("col_mix"     ,2,this->get_nflav(),ncol,nlay); // combination of major species's column amounts
                                                                       // index(1) : reference temperature level
                                                                       // index(2) : flavor
                                                                       // index(3) : layer
    real5d fminor      ("fminor"      ,2,2,this->get_nflav(),ncol,nlay); // interpolation fractions for minor species
                                                                         // index(1) : reference eta level (temperature dependent)
                                                                         // index(2) : reference temperature level
                                                                         // index(3) : flavor
                                                                         // index(4) : layer
    // Error checking
    // Check for initialization
    if (! this->is_initialized()) { stoprun("ERROR: spectral configuration not loaded"); }
    // Check for presence of key species in ty_gas_concs; return error if any key species are not present
    this->check_key_species_present(gas_desc);
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if ( any(play < this->press_ref_min) || any(play > this->press_ref_max) ) {
        stoprun("gas_optics(): array play has values outside range");
      }
      if ( any(plev < this->press_ref_min) || any(plev > this->press_ref_max) ) {
        stoprun("gas_optics(): array plev has values outside range");
      }
      if ( any(tlay < this->temp_ref_min) || any(tlay > this->temp_ref_max) ) {
        stoprun("gas_optics(): array tlay has values outside range");
      }
    #endif
    if (allocated(col_dry)) {
      if (size(col_dry,1) != ncol || size(col_dry,2) != nlay) { stoprun("gas_optics(): array col_dry has wrong size"); }
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (any(col_dry < 0.)) { stoprun("gas_optics(): array col_dry has values outside range"); }
      #endif
    }

    bool use_rayl = allocated(this->krayl);

    int ngas  = this->get_ngas();
    int nflav = this->get_nflav();
    int neta  = this->get_neta();
    int npres = this->get_npres();
    int ntemp = this->get_ntemp();
    // number of minor contributors, total num absorption coeffs
    int nminorlower  = size(this->minor_scales_with_density_lower, 1);
    int nminorklower = size(this->kminor_lower, 1);
    int nminorupper  = size(this->minor_scales_with_density_upper, 1);
    int nminorkupper = size(this->kminor_upper, 1);
    // Fill out the array of volume mixing ratios
    for (int igas = 1 ; igas <= ngas ; igas++) {
      // Get vmr if  gas is provided in ty_gas_concs
      for (size_t igas2 = 0 ; igas2 < gas_desc.gas_name.size() ; igas2++) {
        if ( lower_case(this->gas_names(igas)) == lower_case(gas_desc.gas_name[igas2]) ) {
           real2d vmr_slice = vmr.slice<2>(COLON,COLON,igas);
           gas_desc.get_vmr(this->gas_names(igas), vmr_slice);
        }
      }
    }

    // Compute dry air column amounts (number of molecule per cm^2) if user hasn't provided them
    int idx_h2o = string_loc_in_array("h2o", this->gas_names);
    real2d col_dry_wk;
    if (allocated(col_dry)) {
      col_dry_wk = col_dry;
    } else {
      real2d col_dry_arr = this->get_col_dry(vmr.slice<2>(COLON,COLON,idx_h2o),plev); // dry air column amounts computation
      col_dry_wk = col_dry_arr;
    }
    // compute column gas amounts [molec/cm^2]
    // do ilay = 1, nlay
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
      col_gas(icol,ilay,0) = col_dry_wk(icol,ilay);
    });
    // do igas = 1, ngas
    //   do ilay = 1, nlay
    //     do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<3>(ngas,nlay,ncol) , YAKL_LAMBDA (int igas, int ilay, int icol) {
      col_gas(icol,ilay,igas) = vmr(icol,ilay,igas) * col_dry_wk(icol,ilay);
    });
    // ---- calculate gas optical depths ----
    tau = 0;

    interpolation(ncol, nlay, ngas, nflav, neta, npres, ntemp, this->flavor, this->press_ref_log, this->temp_ref,
                  this->press_ref_log_delta, this->temp_ref_min, this->temp_ref_delta, this->press_ref_trop_log,
                  this->vmr_ref, play, tlay, col_gas, jtemp, fmajor, fminor, col_mix, tropo, jeta, jpress);

    compute_tau_absorption(this->max_gpt_diff_lower, this->max_gpt_diff_upper, ncol, nlay, nband, ngpt, ngas, nflav, neta, npres, ntemp, nminorlower, nminorklower,
                           nminorupper, nminorkupper, idx_h2o, this->gpoint_flavor, this->get_band_lims_gpoint(),
                           this->kmajor, this->kminor_lower, this->kminor_upper, this->minor_limits_gpt_lower,
                           this->minor_limits_gpt_upper, this->minor_scales_with_density_lower,
                           this->minor_scales_with_density_upper, this->scale_by_complement_lower,
                           this->scale_by_complement_upper, this->idx_minor_lower, this->idx_minor_upper,
                           this->idx_minor_scaling_lower, this->idx_minor_scaling_upper, this->kminor_start_lower,
                           this->kminor_start_upper, tropo, col_mix, fmajor, fminor, play, tlay, col_gas,
                           jeta, jtemp, jpress, tau, top_at_1);

    if (allocated(this->krayl)) {
      compute_tau_rayleigh( ncol, nlay, nband, ngpt, ngas, nflav, neta, npres, ntemp, this->gpoint_flavor,
                            this->get_band_lims_gpoint(), this->krayl, idx_h2o, col_dry_wk, col_gas,
                            fminor, jeta, tropo, jtemp, tau_rayleigh);
    }
    combine_and_reorder(tau, tau_rayleigh, allocated(this->krayl), optical_props);
  }



  // Compute Planck source functions at layer centers and levels
  void source(bool top_at_1, int ncol, int nlay, int nbnd, int ngpt, real2d const &play, real2d const &plev, real2d const &tlay,
              real1d const &tsfc, int2d const &jtemp, int2d const &jpress, int4d const &jeta, bool2d const &tropo,
              real6d const &fmajor, SourceFuncLW &sources, real2d const &tlev=real2d()) {
    using yakl::intrinsics::merge;
    using yakl::intrinsics::allocated;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    real3d lay_source_t    ("lay_source_t    ",ngpt,nlay,ncol);
    real3d lev_source_inc_t("lev_source_inc_t",ngpt,nlay,ncol);
    real3d lev_source_dec_t("lev_source_dec_t",ngpt,nlay,ncol);
    real2d sfc_source_t    ("sfc_source_t    ",ngpt     ,ncol);
    // Variables for temperature at layer edges [K] (ncol, nlay+1)
    real2d tlev_arr("tlev_arr",ncol,nlay+1);

    // Source function needs temperature at interfaces/levels and at layer centers
    real2d tlev_wk;
    if (allocated(tlev)) {
      //   Users might have provided these
      tlev_wk = tlev;
    } else {
      tlev_wk = real2d("tlev_wk",ncol,nlay+1);
      // Interpolate temperature to levels if not provided
      //   Interpolation and extrapolation at boundaries is weighted by pressure
      // do ilay = 1, nlay+1
      //   do icol = 1, ncol
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlay+1,ncol) , YAKL_LAMBDA (int ilay, int icol) {
        if (ilay == 1) {
          tlev_wk(icol,1) = tlay(icol,1) + (plev(icol,1)-play(icol,1))*(tlay(icol,2)-tlay(icol,1)) / (play(icol,2)-play(icol,1));
        }
        else if (ilay == nlay+1) {
          tlev_wk(icol,nlay+1) = tlay(icol,nlay) + (plev(icol,nlay+1)-play(icol,nlay))*(tlay(icol,nlay)-tlay(icol,nlay-1)) /
                                                   (play(icol,nlay)-play(icol,nlay-1));
        }
        else {
          tlev_wk(icol,ilay) = ( play(icol,ilay-1)*tlay(icol,ilay-1)*(plev(icol,ilay  )-play(icol,ilay))  +
                                 play(icol,ilay  )*tlay(icol,ilay  )*(play(icol,ilay-1)-plev(icol,ilay)) ) /
                             (plev(icol,ilay)*(play(icol,ilay-1) - play(icol,ilay)));
        }
      });
    }
    // Compute internal (Planck) source functions at layers and levels,
    //  which depend on mapping from spectral space that creates k-distribution.
    int nlayTmp = merge( nlay , 1 , top_at_1 );
    compute_Planck_source(ncol, nlay, nbnd, ngpt, this->get_nflav(), this->get_neta(), this->get_npres(), this->get_ntemp(),
                          this->get_nPlanckTemp(), tlay, tlev_wk, tsfc, nlayTmp, fmajor, jeta, tropo, jtemp, jpress,
                          this->get_gpoint_bands(), this->get_band_lims_gpoint(), this->planck_frac, this->temp_ref_min,
                          this->totplnk_delta, this->totplnk, this->gpoint_flavor, sfc_source_t, lay_source_t, lev_source_inc_t,
                          lev_source_dec_t);

    auto &sources_sfc_source = sources.sfc_source;
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      sources_sfc_source(icol,igpt) = sfc_source_t(igpt,icol);
    });
    reorder123x321(ngpt, nlay, ncol, lay_source_t    , sources.lay_source    );
    reorder123x321(ngpt, nlay, ncol, lev_source_inc_t, sources.lev_source_inc);
    reorder123x321(ngpt, nlay, ncol, lev_source_dec_t, sources.lev_source_dec);
  }



  // Utility function, provided for user convenience
  // computes column amounts of dry air using hydrostatic equation
  real2d get_col_dry(real2d const &vmr_h2o, real2d const &plev, real1d const &latitude=real1d()) {
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;
    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;

    // first and second term of Helmert formula
    real constexpr helmert1 = 9.80665;
    real constexpr helmert2 = 0.02586;
    int ncol = size(plev,1);
    int nlev = size(plev,2);
    real1d g0("g0",size(plev,1));
    if (allocated(latitude)) {
      // A purely OpenACC implementation would probably compute g0 within the kernel below
      // do icol = 1, ncol
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<1>(ncol) , YAKL_LAMBDA (int icol) {
        g0(icol) = helmert1 - helmert2 * cos(2.0 * M_PI * latitude(icol) / 180.0); // acceleration due to gravity [m/s^2]
      });
    } else {
      // do icol = 1, ncol
      YAKL_SCOPE( grav , ::grav );
      parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<1>(ncol) , YAKL_LAMBDA (int icol) {
        g0(icol) = grav;
      });
    }

    real2d col_dry("col_dry",size(plev,1),size(plev,2)-1);
    // do ilev = 1, nlev-1
    //   do icol = 1, ncol
    YAKL_SCOPE( m_dry , ::m_dry );
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<2>(nlev-1,ncol) , YAKL_LAMBDA (int ilev , int icol) {
      real delta_plev = abs(plev(icol,ilev) - plev(icol,ilev+1));
      // Get average mass of moist air per mole of moist air
      real fact = 1. / (1.+vmr_h2o(icol,ilev));
      real m_air = (m_dry + m_h2o * vmr_h2o(icol,ilev)) * fact;
      col_dry(icol,ilev) = 10. * delta_plev * avogad * fact/(1000.*m_air*100.*g0(icol));
    });
    return col_dry;
  }



 // Utility function to combine optical depths from gas absorption and Rayleigh scattering
 //   (and reorder them for convenience, while we're at it)
 void combine_and_reorder(real3d const &tau, real3d const &tau_rayleigh, bool has_rayleigh, OpticalProps1scl &optical_props) {
    using yakl::intrinsics::size;

    int ncol = size(tau,3);
    int nlay = size(tau,2);
    int ngpt = size(tau,1);
    reorder123x321(ngpt, nlay, ncol, tau, optical_props.tau);
  }



 // Utility function to combine optical depths from gas absorption and Rayleigh scattering
 //   (and reorder them for convenience, while we're at it)
 void combine_and_reorder(real3d const &tau, real3d const &tau_rayleigh, bool has_rayleigh, OpticalProps2str &optical_props) {
    using yakl::intrinsics::size;

    int ncol = size(tau,3);
    int nlay = size(tau,2);
    int ngpt = size(tau,1);
    if (has_rayleigh) {
      // combine optical depth and rayleigh scattering
      combine_and_reorder_2str(ncol, nlay, ngpt, tau, tau_rayleigh, optical_props.tau, optical_props.ssa, optical_props.g);
    } else {
      // index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
      reorder123x321(ngpt, nlay, ncol, tau, optical_props.tau);
      optical_props.ssa = 0;
      optical_props.g   = 0;
    }
  }



  void print_norms() const {
    using yakl::intrinsics::count;
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

                                                      std::cout << "name                                  : " << std::setw(20) << name                                   << "\n";
                                                      std::cout << "totplnk_delta                         : " << std::setw(20) << totplnk_delta                          << "\n";
                                                      std::cout << "press_ref_min                         : " << std::setw(20) << press_ref_min                          << "\n";
                                                      std::cout << "press_ref_max                         : " << std::setw(20) << press_ref_max                          << "\n";
                                                      std::cout << "temp_ref_min                          : " << std::setw(20) << temp_ref_min                           << "\n";
                                                      std::cout << "temp_ref_max                          : " << std::setw(20) << temp_ref_max                           << "\n";
                                                      std::cout << "press_ref_log_delta                   : " << std::setw(20) << press_ref_log_delta                    << "\n";
                                                      std::cout << "temp_ref_delta                        : " << std::setw(20) << temp_ref_delta                         << "\n";
                                                      std::cout << "press_ref_trop_log                    : " << std::setw(20) << press_ref_trop_log                     << "\n";
    if (allocated(gas_names                      )) { std::cout << "gas_names                             : " << std::setw(20) << gas_names                              << "\n"; }
    if (allocated(band2gpt                       )) { std::cout << "sum(band2gpt     )                    : " << std::setw(20) << sum(band2gpt     )                     << "\n"; }
    if (allocated(gpt2band                       )) { std::cout << "sum(gpt2band     )                    : " << std::setw(20) << sum(gpt2band     )                     << "\n"; }
    if (allocated(band_lims_wvn                  )) { std::cout << "sum(band_lims_wvn)                    : " << std::setw(20) << sum(band_lims_wvn)                     << "\n"; }
    if (allocated(press_ref                      )) { std::cout << "sum(press_ref    )                    : " << std::setw(20) << sum(press_ref    )                     << "\n"; }
    if (allocated(press_ref_log                  )) { std::cout << "sum(press_ref_log)                    : " << std::setw(20) << sum(press_ref_log)                     << "\n"; }
    if (allocated(temp_ref                       )) { std::cout << "sum(temp_ref     )                    : " << std::setw(20) << sum(temp_ref     )                     << "\n"; }
    if (allocated(vmr_ref                        )) { std::cout << "sum(vmr_ref                )          : " << std::setw(20) << sum(vmr_ref                )           << "\n"; }
    if (allocated(flavor                         )) { std::cout << "sum(flavor                 )          : " << std::setw(20) << sum(flavor                 )           << "\n"; }
    if (allocated(gpoint_flavor                  )) { std::cout << "sum(gpoint_flavor          )          : " << std::setw(20) << sum(gpoint_flavor          )           << "\n"; }
    if (allocated(kmajor                         )) { std::cout << "sum(kmajor                 )          : " << std::setw(20) << sum(kmajor                 )           << "\n"; }
    if (allocated(minor_limits_gpt_lower         )) { std::cout << "sum(minor_limits_gpt_lower )          : " << std::setw(20) << sum(minor_limits_gpt_lower )           << "\n"; }
    if (allocated(minor_limits_gpt_upper         )) { std::cout << "sum(minor_limits_gpt_upper )          : " << std::setw(20) << sum(minor_limits_gpt_upper )           << "\n"; }
    if (allocated(idx_minor_lower                )) { std::cout << "sum(idx_minor_lower        )          : " << std::setw(20) << sum(idx_minor_lower        )           << "\n"; }
    if (allocated(idx_minor_upper                )) { std::cout << "sum(idx_minor_upper        )          : " << std::setw(20) << sum(idx_minor_upper        )           << "\n"; }
    if (allocated(idx_minor_scaling_lower        )) { std::cout << "sum(idx_minor_scaling_lower)          : " << std::setw(20) << sum(idx_minor_scaling_lower)           << "\n"; }
    if (allocated(idx_minor_scaling_upper        )) { std::cout << "sum(idx_minor_scaling_upper)          : " << std::setw(20) << sum(idx_minor_scaling_upper)           << "\n"; }
    if (allocated(kminor_start_lower             )) { std::cout << "sum(kminor_start_lower     )          : " << std::setw(20) << sum(kminor_start_lower     )           << "\n"; }
    if (allocated(kminor_start_upper             )) { std::cout << "sum(kminor_start_upper     )          : " << std::setw(20) << sum(kminor_start_upper     )           << "\n"; }
    if (allocated(kminor_lower                   )) { std::cout << "sum(kminor_lower           )          : " << std::setw(20) << sum(kminor_lower           )           << "\n"; }
    if (allocated(kminor_upper                   )) { std::cout << "sum(kminor_upper           )          : " << std::setw(20) << sum(kminor_upper           )           << "\n"; }
    if (allocated(krayl                          )) { std::cout << "sum(krayl                  )          : " << std::setw(20) << sum(krayl                  )           << "\n"; }
    if (allocated(planck_frac                    )) { std::cout << "sum(planck_frac            )          : " << std::setw(20) << sum(planck_frac            )           << "\n"; }
    if (allocated(totplnk                        )) { std::cout << "sum(totplnk                )          : " << std::setw(20) << sum(totplnk                )           << "\n"; }
    if (allocated(solar_src                      )) { std::cout << "sum(solar_src              )          : " << std::setw(20) << sum(solar_src              )           << "\n"; }
    if (allocated(minor_scales_with_density_lower)) { std::cout << "count(minor_scales_with_density_lower): " << std::setw(20) << count(minor_scales_with_density_lower) << "\n"; }
    if (allocated(minor_scales_with_density_upper)) { std::cout << "count(minor_scales_with_density_upper): " << std::setw(20) << count(minor_scales_with_density_upper) << "\n"; }
    if (allocated(scale_by_complement_lower      )) { std::cout << "count(scale_by_complement_lower      ): " << std::setw(20) << count(scale_by_complement_lower      ) << "\n"; }
    if (allocated(scale_by_complement_upper      )) { std::cout << "count(scale_by_complement_upper      ): " << std::setw(20) << count(scale_by_complement_upper      ) << "\n"; }
    if (allocated(is_key                         )) { std::cout << "count(is_key                         ): " << std::setw(20) << count(is_key                         ) << "\n"; }
  }


};
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
class GasOpticsRRTMGPK : public OpticalPropsK {
public:
  // RRTMGP computes absorption in each band arising from
  //   two major species in each band, which are combined to make
  //     a relative mixing ratio eta and a total column amount (col_mix)
  //   contributions from zero or more minor species whose concentrations
  //     may be scaled by other components of the atmosphere
  // Absorption coefficients are interpolated from tables on a pressure/temperature/(eta) grid
  // Interpolation variables: Temperature and pressure grids
  real1dk press_ref;
  real1dk press_ref_log;
  real1dk temp_ref;

  // Derived and stored for convenience:
  //   Min and max for temperature and pressure intepolation grids
  //   difference in ln pressure between consecutive reference levels
  //   log of reference pressure separating the lower and upper atmosphere
  real press_ref_min;
  real press_ref_max;
  real temp_ref_min;
  real temp_ref_max;
  real press_ref_log_delta;
  real temp_ref_delta;
  real press_ref_trop_log;

  int max_gpt_diff_lower;
  int max_gpt_diff_upper;

  // Major absorbers ("key species")
  //   Each unique set of major species is called a flavor.
  // Names  and reference volume mixing ratios of major gases
  string1dv gas_names; // gas names
  realOff3dk vmr_ref;  // vmr_ref(lower or upper atmosphere, gas, temp)

  // Which two gases are in each flavor? By index
  int2dk flavor;        // major species pair; (2,nflav)

  // Which flavor for each g-point? One each for lower, upper atmosphere
  int2dk gpoint_flavor; // flavor = gpoint_flavor(2, g-point)

  // Major gas absorption coefficients
  real4dk kmajor;       //  kmajor(g-point,eta,pressure,temperature)

  // Minor species, independently for upper and lower atmospheres
  //   Array extents in the n_minor dimension will differ between upper and lower atmospheres
  //   Each contribution has starting and ending g-points
  int2dk minor_limits_gpt_lower;
  int2dk minor_limits_gpt_upper;

  // Minor gas contributions might be scaled by other gas amounts; if so we need to know
  //   the total density and whether the contribution is scaled by the partner gas
  //   or its complement (i.e. all other gases)
  // Water vapor self- and foreign continua work like this, as do
  //   all collision-induced abosption pairs
  bool1dk minor_scales_with_density_lower;
  bool1dk minor_scales_with_density_upper;
  bool1dk scale_by_complement_lower;
  bool1dk scale_by_complement_upper;
  int1dk idx_minor_lower;
  int1dk idx_minor_upper;
  int1dk idx_minor_scaling_lower;
  int1dk idx_minor_scaling_upper;

  // Index into table of absorption coefficients
  int1dk kminor_start_lower;
  int1dk kminor_start_upper;

  // The absorption coefficients themselves
  real3dk kminor_lower; // kminor_lower(n_minor,eta,temperature)
  real3dk kminor_upper; // kminor_upper(n_minor,eta,temperature)

  // Rayleigh scattering coefficients
  real4dk krayl; // krayl(g-point,eta,temperature,upper/lower atmosphere)

  // Planck function spectral mapping
  //   Allocated only when gas optics object is internal-source
  real4dk planck_frac;   // stored fraction of Planck irradiance in band for given g-point
                        // planck_frac(g-point, eta, pressure, temperature)
  real2dk totplnk;       // integrated Planck irradiance by band; (Planck temperatures,band)
  real   totplnk_delta; // temperature steps in totplnk

  // Solar source function spectral mapping
  //   Allocated only when gas optics object is external-source
  real1dk solar_src; // incoming solar irradiance(g-point)

  // Ancillary
  // Index into %gas_names -- is this a key species in any band?
  bool1dk is_key;

  void finalize () {
      press_ref = decltype(press_ref)();
      press_ref_log = decltype(press_ref_log)();
      temp_ref = decltype(temp_ref)();
      gas_names = decltype(gas_names)();  // gas names
      vmr_ref = decltype(vmr_ref)();      // vmr_ref(lower or upper atmosphere, gas, temp)
      flavor = decltype(flavor)();        // major species pair; (2,nflav)
      gpoint_flavor = decltype(gpoint_flavor)(); // flavor = gpoint_flavor(2, g-point)
      kmajor = decltype(kmajor)();       //  kmajor(g-point,eta,pressure,temperature)
      minor_limits_gpt_lower = decltype(minor_limits_gpt_lower)();
      minor_limits_gpt_upper = decltype(minor_limits_gpt_upper)();
      minor_scales_with_density_lower = decltype(minor_scales_with_density_lower)();
      minor_scales_with_density_upper = decltype(minor_scales_with_density_upper)();
      scale_by_complement_lower = decltype(scale_by_complement_lower)();
      scale_by_complement_upper = decltype(scale_by_complement_upper)();
      idx_minor_lower = decltype(idx_minor_lower)();
      idx_minor_upper = decltype(idx_minor_upper)();
      idx_minor_scaling_lower = decltype(idx_minor_scaling_lower)();
      idx_minor_scaling_upper = decltype(idx_minor_scaling_upper)();
      kminor_start_lower = decltype(kminor_start_lower)();
      kminor_start_upper = decltype(kminor_start_upper)();
      kminor_lower = decltype(kminor_lower)(); // kminor_lower(n_minor,eta,temperature)
      kminor_upper = decltype(kminor_upper)(); // kminor_upper(n_minor,eta,temperature)
      krayl = decltype(krayl)(); // krayl(g-point,eta,temperature,upper/lower atmosphere)
      planck_frac = decltype(planck_frac)();   // stored fraction of Planck irradiance in band for given g-point
      totplnk = decltype(totplnk)();       // integrated Planck irradiance by band; (Planck temperatures,band)
      solar_src = decltype(solar_src)(); // incoming solar irradiance(g-point)
      is_key = decltype(is_key)();

      // Free memory in base class
      OpticalPropsK::finalize();
  }

  // Everything except GasConcs is on the host by default, and the available_gases.gas_name is on the host as well
  // Things will be copied to the GPU outside of this routine and stored into class variables
  void reduce_minor_arrays(GasConcsK   const &available_gases,
                           string1dv   const &gas_names,
                           string1dv   const &gas_minor,
                           string1dv   const &identifier_minor,
                           realHost3dk const &kminor_atm,
                           string1dv   const &minor_gases_atm,
                           intHost2dk  const &minor_limits_gpt_atm,
                           boolHost1dk const &minor_scales_with_density_atm,
                           string1dv   const &scaling_gas_atm,
                           boolHost1dk const &scale_by_complement_atm,
                           intHost1dk  const &kminor_start_atm,
                           realHost3dk       &kminor_atm_red,
                           string1dv         &minor_gases_atm_red,
                           intHost2dk        &minor_limits_gpt_atm_red,
                           boolHost1dk       &minor_scales_with_density_atm_red,
                           string1dv         &scaling_gas_atm_red,
                           boolHost1dk       &scale_by_complement_atm_red,
                           intHost1dk        &kminor_start_atm_red) {
    int nm = minor_gases_atm.size();  // Size of the larger list of minor gases
    int tot_g = 0;
    int red_nm = 0;                    // Reduced number of minor gasses (only the ones we need)
    boolHost1dk gas_is_present("gas_is_present",nm);   // Determines whether a gas in the list is needed
    // Determine the gasses needed
    for (int i=0; i < nm; i++) {
      int idx_mnr = string_loc_in_array(minor_gases_atm[i], identifier_minor);
      gas_is_present(i) = string_in_array(gas_minor[idx_mnr], available_gases.gas_name);
      if (gas_is_present(i)) {
        tot_g = tot_g + (minor_limits_gpt_atm(1,i)-minor_limits_gpt_atm(0,i)+1);
        red_nm = red_nm + 1;
      }
    }

    // Allocate reduced arrays
    minor_gases_atm_red               = string1dv  (red_nm);
    minor_scales_with_density_atm_red = boolHost1dk("minor_scales_with_density_atm_red"  ,red_nm);
    scaling_gas_atm_red               = string1dv  (red_nm);
    scale_by_complement_atm_red       = boolHost1dk("scale_by_complement_atm_red      "  ,red_nm);
    kminor_start_atm_red              = intHost1dk ("kminor_start_atm_red             "  ,red_nm);
    minor_limits_gpt_atm_red          = intHost2dk ("minor_limits_gpt_atm_red         ",2,red_nm);
    kminor_atm_red                    = realHost3dk("kminor_atm_red                   ",tot_g , kminor_atm.extent(1), kminor_atm.extent(2));

    if (red_nm == nm) {
      // If the gasses listed exactly matches the gasses needed, just copy it
      minor_gases_atm_red = minor_gases_atm;
      scaling_gas_atm_red = scaling_gas_atm;
      Kokkos::deep_copy(kminor_atm_red, kminor_atm);
      Kokkos::deep_copy(minor_limits_gpt_atm_red, minor_limits_gpt_atm);
      Kokkos::deep_copy(minor_scales_with_density_atm_red, minor_scales_with_density_atm);
      Kokkos::deep_copy(scale_by_complement_atm_red, scale_by_complement_atm);
      Kokkos::deep_copy(kminor_start_atm_red, kminor_start_atm);
    } else {
      // Otherwise, pack into reduced arrays
      int slot = 0;
      for (int i=0; i < nm; i++) {
        if (gas_is_present(i)) {
          minor_gases_atm_red              [slot] = minor_gases_atm              [i];
          scaling_gas_atm_red              [slot] = scaling_gas_atm              [i];
          minor_scales_with_density_atm_red(slot) = minor_scales_with_density_atm(i);
          scale_by_complement_atm_red      (slot) = scale_by_complement_atm      (i);
          kminor_start_atm_red             (slot) = kminor_start_atm             (i);
          slot++;
        }
      }

      slot = -1;
      int n_elim = 0;
      for (int i=0 ; i < nm ; i++) {
        int ng = minor_limits_gpt_atm(1,i)-minor_limits_gpt_atm(0,i)+1;
        if (gas_is_present(i)) {
          ++slot;
          minor_limits_gpt_atm_red   (0,slot) = minor_limits_gpt_atm(0,i);
          minor_limits_gpt_atm_red   (1,slot) = minor_limits_gpt_atm(1,i);
          kminor_start_atm_red         (slot) = kminor_start_atm(i)-n_elim;
          for (int l=0 ; l < kminor_atm.extent(2); l++ ) {
            for (int k=0 ; k < kminor_atm.extent(1) ; k++ ) {
              for (int j=0 ; j < ng ; j++) {
                kminor_atm_red(kminor_start_atm_red(slot)+j,k,l) = kminor_atm(kminor_start_atm(i)+j,k,l);
              }
            }
          }
        } else {
          n_elim = n_elim + ng;
        }
      }
    }
  }

  // create index list for extracting col_gas needed for minor gas optical depth calculations
  void create_idx_minor(string1dv const &gas_names, string1dv const &gas_minor, string1dv const &identifier_minor,
                        string1dv const &minor_gases_atm, intHost1dk &idx_minor_atm) {
    idx_minor_atm = intHost1dk("idx_minor_atm", minor_gases_atm.size());
    for (size_t imnr=0 ; imnr < minor_gases_atm.size() ; imnr++) {
      // Find identifying string for minor species in list of possible identifiers (e.g. h2o_slf)
      int idx_mnr     = string_loc_in_array(minor_gases_atm[imnr], identifier_minor);
      // Find name of gas associated with minor species identifier (e.g. h2o)
      idx_minor_atm(imnr) = string_loc_in_array(gas_minor[idx_mnr], gas_names);
    }
  }

  // create index for special treatment in density scaling of minor gases
  void create_idx_minor_scaling(string1dv const &gas_names, string1dv const &scaling_gas_atm,
                                intHost1dk &idx_minor_scaling_atm) {
    idx_minor_scaling_atm = intHost1dk("idx_minor_scaling_atm", scaling_gas_atm.size());
    for (auto imnr=0 ; imnr < scaling_gas_atm.size() ; imnr++) {
      // This will be -1 if there's no interacting gas
      idx_minor_scaling_atm(imnr) = string_loc_in_array(scaling_gas_atm[imnr], gas_names);
    }
  }

  void create_key_species_reduce(string1dv const &gas_names, string1dv const &gas_names_red, intHost3dk const &key_species,
                                 intHost3dk &key_species_red, boolHost1dk &key_species_present_init) {
    int np = key_species.extent(0);
    int na = key_species.extent(1);
    int nt = key_species.extent(2);
    key_species_red = intHost3dk("key_species_red",np,na,nt);
    key_species_present_init = boolHost1dk("key_species_present_init", gas_names.size());
    Kokkos::deep_copy(key_species_present_init, true);

    for (int ip=0 ; ip < np ; ip++) {
      for (int ia=0 ; ia < na ; ia++) {
        for (int it=0 ; it < nt ; it++) {
          if (key_species(ip,ia,it) != -1) {
            key_species_red(ip,ia,it) = string_loc_in_array(gas_names[key_species(ip,ia,it)],gas_names_red);
            if (key_species_red(ip,ia,it) == -1) {
              key_species_present_init(key_species(ip,ia,it)) = false;
            }
          } else {
            key_species_red(ip,ia,it) = -1;
          }
        }
      }
    }
  }

  // Create flavor list
  // An unordered array of extent (2,:) containing all possible pairs of key species used in either upper or lower atmos
  void create_flavor(intHost3dk const &key_species, intHost2dk &flavor) {
    // prepare list of key_species
    int i = 0;
    intHost2dk key_species_list("key_species_list", 2, key_species.extent(2)*2);
    for (int ibnd=0 ; ibnd < key_species.extent(2) ; ibnd++) {
      for (int iatm=0 ; iatm < key_species.extent(0) ; iatm++, i++) {
        key_species_list(0,i) = key_species(0,iatm,ibnd);
        key_species_list(1,i) = key_species(1,iatm,ibnd);
      }
    }
    // rewrite single key_species pairs
    for (int i=0 ; i < key_species_list.extent(1) ; i++) {
      if (key_species_list(0,i) == -1 && key_species_list(1,i) == -1) {
        key_species_list(0,i) = 1;
        key_species_list(1,i) = 1;
      }
    }
    // count unique key species pairs
    int iflavor = 0;
    for (int i=0; i < key_species_list.extent(1) ; i++) {
      // Loop through previous pairs. Only increment iflavor if we haven't seen this pair before
      bool unique = true;
      for (int j=0; j <= i-1 ; j++) {
        if ( key_species_list(0,j) == key_species_list(0,i) && key_species_list(1,j) == key_species_list(1,i) ) {
          unique = false;
          break;
        }
      }
      if (unique) { ++iflavor; }
    }
    // fill flavors
    flavor = intHost2dk("flavor",2,iflavor);
    iflavor = 0;
    for (int i=0 ; i < key_species_list.extent(1) ; i++) {
      bool unique = true;
      for (int j=0; j <= i-1 ; j++) {
        if ( key_species_list(0,j) == key_species_list(0,i) && key_species_list(1,j) == key_species_list(1,i) ) {
          unique = false;
          break;
        }
      }
      if (unique) {
        flavor(0,iflavor) = key_species_list(0,i);
        flavor(1,iflavor) = key_species_list(1,i);
        ++iflavor;
      }
    }
  }

  // create gpoint_flavor list: a map pointing from each g-point to the corresponding entry in the "flavor list"
  void create_gpoint_flavor(intHost3dk const &key_species, intHost1dk const &gpt2band, intHost2dk const &flavor,
                            intHost2dk &gpoint_flavor) {
    int ngpt = gpt2band.extent(0);
    gpoint_flavor = intHost2dk("gpoint_flavor",2,ngpt);
    for (int igpt=0 ; igpt < ngpt ; igpt++) {
      for (int iatm = 0 ; iatm < 2 ; iatm++) {
        int key_species_pair2flavor = -1;
        for (int iflav=0 ; iflav < flavor.extent(1) ; iflav++) {
          int key_species1 = key_species(0,iatm,gpt2band(igpt));
          int key_species2 = key_species(1,iatm,gpt2band(igpt));
          if (key_species1 == -1 && key_species2 == -1) {
            key_species1 = 1;
            key_species2 = 1;
          }
          if ( flavor(0,iflav) == key_species1 && flavor(1,iflav) == key_species2 ) {
            key_species_pair2flavor = iflav;
          }
        }
        gpoint_flavor(iatm,igpt) = key_species_pair2flavor;
      }
    }
  }

  // Initialize absorption coefficient arrays,
  //   including Rayleigh scattering tables if provided (allocated)
  void init_abs_coeffs(GasConcsK   const &available_gases,
                       string1dv   const &gas_names,
                       intHost3dk  const &key_species,
                       intHost2dk  const &band2gpt,
                       realHost2dk const &band_lims_wavenum,
                       realHost1dk const &press_ref,
                       realHost1dk const &temp_ref,
                       real              press_ref_trop,
                       real              temp_ref_p,
                       real              temp_ref_t,
                       realHost3dk const &vmr_ref,
                       realHost4dk const &kmajor,
                       realHost3dk const &kminor_lower,
                       realHost3dk const &kminor_upper,
                       string1dv   const &gas_minor,
                       string1dv   const &identifier_minor,
                       string1dv   const &minor_gases_lower,
                       string1dv   const &minor_gases_upper,
                       intHost2dk  const &minor_limits_gpt_lower,
                       intHost2dk  const &minor_limits_gpt_upper,
                       boolHost1dk const &minor_scales_with_density_lower,
                       boolHost1dk const &minor_scales_with_density_upper,
                       string1dv   const &scaling_gas_lower,
                       string1dv   const &scaling_gas_upper,
                       boolHost1dk const &scale_by_complement_lower,
                       boolHost1dk const &scale_by_complement_upper,
                       intHost1dk  const &kminor_start_lower,
                       intHost1dk  const &kminor_start_upper,
                       realHost3dk const &rayl_lower,
                       realHost3dk const &rayl_upper) {
    auto band_lims_wavenum_h = Kokkos::create_mirror_view(DefaultDevice(), band_lims_wavenum);
    auto band2gpt_h = Kokkos::create_mirror_view(DefaultDevice(), band2gpt);
    Kokkos::deep_copy(band_lims_wavenum_h, band_lims_wavenum);
    Kokkos::deep_copy(band2gpt_h, band2gpt);
    OpticalPropsK::init(band_lims_wavenum_h, band2gpt_h);

    // Which gases known to the gas optics are present in the host model (available_gases)?
    for (auto item : gas_names) {
      if (string_in_array(item, available_gases.gas_name)) {
        this->gas_names.push_back(item);
      }
    }
    // Now the number of gases is the union of those known to the k-distribution and provided by the host model

    // Initialize the gas optics object, keeping only those gases known to the gas optics and also present in the host model
    const int ngas = this->gas_names.size();
    const int vmr_e0 = vmr_ref.extent(0);
    const int vmr_e2 = vmr_ref.extent(2);

    realOffHost3dk vmr_ref_red("vmr_ref_red", std::make_pair(0, vmr_e0-1), std::make_pair(-1, ngas-1), std::make_pair(0, vmr_e2-1));
    // Gas 0 is used in single-key species method, set to 1.0 (col_dry)
    for (int k=0 ; k < vmr_e2 ; k++) {
      for (int j=0 ; j < vmr_e0 ; j++) {
        vmr_ref_red(j,-1,k) = vmr_ref(j,0,k);
      }
    }
    for (int i=0 ; i < ngas ; i++) {
      int idx = string_loc_in_array(this->gas_names[i], gas_names);
      for (int k=0 ; k < vmr_e2 ; k++) {
        for (int j=0 ; j < vmr_e0 ; j++) {
          vmr_ref_red(j,i,k) = vmr_ref(j,idx+1,k);
        }
      }
    }
    // Allocate class copy, and deep copy to the class data member
    this->vmr_ref = Kokkos::create_mirror_view(DefaultDevice(), vmr_ref_red);
    Kokkos::deep_copy(this->vmr_ref, vmr_ref_red);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // REDUCE MINOR ARRAYS SO VARIABLES ONLY CONTAIN MINOR GASES THAT ARE AVAILABLE
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // LOWER MINOR GASSES
    string1dv   minor_gases_lower_red;
    string1dv   scaling_gas_lower_red;
    realHost3dk kminor_lower_red;
    intHost2dk  minor_limits_gpt_lower_red;
    boolHost1dk minor_scales_with_density_lower_red;
    boolHost1dk scale_by_complement_lower_red;
    intHost1dk  kminor_start_lower_red;

    reduce_minor_arrays(available_gases, gas_names, gas_minor, identifier_minor, kminor_lower, minor_gases_lower,
                        minor_limits_gpt_lower, minor_scales_with_density_lower, scaling_gas_lower,
                        scale_by_complement_lower, kminor_start_lower, kminor_lower_red, minor_gases_lower_red,
                        minor_limits_gpt_lower_red, minor_scales_with_density_lower_red, scaling_gas_lower_red,
                        scale_by_complement_lower_red, kminor_start_lower_red);

    this->kminor_lower                    = Kokkos::create_mirror_view(DefaultDevice(), kminor_lower_red);
    this->minor_limits_gpt_lower          = Kokkos::create_mirror_view(DefaultDevice(), minor_limits_gpt_lower_red);
    this->minor_scales_with_density_lower = Kokkos::create_mirror_view(DefaultDevice(), minor_scales_with_density_lower_red);
    this->scale_by_complement_lower       = Kokkos::create_mirror_view(DefaultDevice(), scale_by_complement_lower_red);
    this->kminor_start_lower              = Kokkos::create_mirror_view(DefaultDevice(), kminor_start_lower_red);

    Kokkos::deep_copy(this->kminor_lower, kminor_lower_red);
    Kokkos::deep_copy(this->minor_limits_gpt_lower, minor_limits_gpt_lower_red);
    Kokkos::deep_copy(this->minor_scales_with_density_lower, minor_scales_with_density_lower_red);
    Kokkos::deep_copy(this->scale_by_complement_lower, scale_by_complement_lower_red);
    Kokkos::deep_copy(this->kminor_start_lower, kminor_start_lower_red);

    // Find the largest number of g-points per band
    this->max_gpt_diff_lower = std::numeric_limits<int>::lowest();
    for (int i=0; i<minor_limits_gpt_lower_red.extent(1); i++) {
      this->max_gpt_diff_lower = std::max( this->max_gpt_diff_lower , minor_limits_gpt_lower_red(1,i) - minor_limits_gpt_lower_red(0,i) );
    }

    // UPPER MINOR GASSES
    string1dv minor_gases_upper_red;
    string1dv scaling_gas_upper_red;
    realHost3dk kminor_upper_red;
    intHost2dk  minor_limits_gpt_upper_red;
    boolHost1dk minor_scales_with_density_upper_red;
    boolHost1dk scale_by_complement_upper_red;
    intHost1dk  kminor_start_upper_red;

    reduce_minor_arrays(available_gases, gas_names, gas_minor, identifier_minor, kminor_upper, minor_gases_upper,
                        minor_limits_gpt_upper, minor_scales_with_density_upper, scaling_gas_upper,
                        scale_by_complement_upper, kminor_start_upper, kminor_upper_red, minor_gases_upper_red,
                        minor_limits_gpt_upper_red, minor_scales_with_density_upper_red, scaling_gas_upper_red,
                        scale_by_complement_upper_red, kminor_start_upper_red);

    this->kminor_upper                    = Kokkos::create_mirror_view(DefaultDevice(), kminor_upper_red);
    this->minor_limits_gpt_upper          = Kokkos::create_mirror_view(DefaultDevice(), minor_limits_gpt_upper_red);
    this->minor_scales_with_density_upper = Kokkos::create_mirror_view(DefaultDevice(), minor_scales_with_density_upper_red);
    this->scale_by_complement_upper       = Kokkos::create_mirror_view(DefaultDevice(), scale_by_complement_upper_red);
    this->kminor_start_upper              = Kokkos::create_mirror_view(DefaultDevice(), kminor_start_upper_red);

    Kokkos::deep_copy(this->kminor_upper, kminor_upper_red);
    Kokkos::deep_copy(this->minor_limits_gpt_upper, minor_limits_gpt_upper_red);
    Kokkos::deep_copy(this->minor_scales_with_density_upper, minor_scales_with_density_upper_red);
    Kokkos::deep_copy(this->scale_by_complement_upper, scale_by_complement_upper_red);
    Kokkos::deep_copy(this->kminor_start_upper, kminor_start_upper_red);

    // Find the largest number of g-points per band
    this->max_gpt_diff_upper = std::numeric_limits<int>::lowest();
    for (int i=0; i<minor_limits_gpt_upper_red.extent(1); i++) {
      this->max_gpt_diff_upper = std::max( this->max_gpt_diff_upper , minor_limits_gpt_upper_red(1,i) - minor_limits_gpt_upper_red(0,i) );
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // HANDLE ARRAYS NOT REDUCED BY THE PRESENCE, OR LACK THEREOF, OF A GAS
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    this->press_ref = Kokkos::create_mirror_view(DefaultDevice(), press_ref);
    this->temp_ref  = Kokkos::create_mirror_view(DefaultDevice(), temp_ref);
    this->kmajor    = Kokkos::create_mirror_view(DefaultDevice(), kmajor);
    Kokkos::deep_copy(this->press_ref, press_ref);
    Kokkos::deep_copy(this->temp_ref, temp_ref);
    Kokkos::deep_copy(this->kmajor, kmajor);

    // Process rayl_lower and rayl_upper into a combined this->krayl
    if (rayl_lower.is_allocated() != rayl_upper.is_allocated()) {
      stoprun("rayl_lower and rayl_upper must have the same allocation status");
    }
    if (rayl_lower.is_allocated()) {
      realHost4dk krayltmp("krayltmp",rayl_lower.extent(0),rayl_lower.extent(1),rayl_lower.extent(2),2);
      for (int k=0 ; k < rayl_lower.extent(2) ; k++ ) {
        for (int j=0 ; j < rayl_lower.extent(1) ; j++ ) {
          for (int i=0 ; i < rayl_lower.extent(0) ; i++ ) {
            krayltmp(i,j,k,0) = rayl_lower(i,j,k);
            krayltmp(i,j,k,1) = rayl_upper(i,j,k);
          }
        }
      }
      this->krayl = Kokkos::create_mirror_view(DefaultDevice(), krayltmp);
      Kokkos::deep_copy(this->krayl, krayltmp);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // POST PROCESSING
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // creates log reference pressure
    this->press_ref_log = real1dk("press_ref_log", this->press_ref.extent(0));
    // Running a kernel because it's more convenient in this case
    Kokkos::parallel_for( this->press_ref.extent(0) , KOKKOS_LAMBDA (int i) {
      this->press_ref_log(i) = log(this->press_ref(i));
    });

    // log scale of reference pressure (this is a scalar, not an array)
    this->press_ref_trop_log = log(press_ref_trop);

    // Get index of gas (if present) for determining col_gas
    intHost1dk idx_minor_lower_tmp;
    intHost1dk idx_minor_upper_tmp;
    create_idx_minor(this->gas_names, gas_minor, identifier_minor, minor_gases_lower_red, idx_minor_lower_tmp);
    create_idx_minor(this->gas_names, gas_minor, identifier_minor, minor_gases_upper_red, idx_minor_upper_tmp);
    this->idx_minor_lower = Kokkos::create_mirror_view(DefaultDevice(), idx_minor_lower_tmp);
    this->idx_minor_upper = Kokkos::create_mirror_view(DefaultDevice(), idx_minor_upper_tmp);
    Kokkos::deep_copy(this->idx_minor_lower, idx_minor_lower_tmp);
    Kokkos::deep_copy(this->idx_minor_upper, idx_minor_upper_tmp);
    // Get index of gas (if present) that has special treatment in density scaling
    intHost1dk idx_minor_scaling_lower_tmp;
    intHost1dk idx_minor_scaling_upper_tmp;
    create_idx_minor_scaling(this->gas_names, scaling_gas_lower_red, idx_minor_scaling_lower_tmp);
    create_idx_minor_scaling(this->gas_names, scaling_gas_upper_red, idx_minor_scaling_upper_tmp);
    this->idx_minor_scaling_lower = Kokkos::create_mirror_view(DefaultDevice(), idx_minor_scaling_lower_tmp);
    this->idx_minor_scaling_upper = Kokkos::create_mirror_view(DefaultDevice(), idx_minor_scaling_upper_tmp);
    Kokkos::deep_copy(this->idx_minor_scaling_lower, idx_minor_scaling_lower_tmp);
    Kokkos::deep_copy(this->idx_minor_scaling_upper, idx_minor_scaling_upper_tmp);

    // create flavor list
    // Reduce (remap) key_species list; checks that all key gases are present in incoming
    boolHost1dk key_species_present_init;
    intHost3dk  key_species_red;
    create_key_species_reduce(gas_names, this->gas_names, key_species, key_species_red, key_species_present_init);
    // create flavor and gpoint_flavor lists
    intHost2dk flavor_tmp;
    intHost2dk gpoint_flavor_tmp;
    auto gpoint_bands_tmp = Kokkos::create_mirror_view(this->get_gpoint_bands());
    Kokkos::deep_copy(gpoint_bands_tmp, this->get_gpoint_bands());
    create_flavor       (key_species_red, flavor_tmp);
    create_gpoint_flavor(key_species_red, gpoint_bands_tmp, flavor_tmp, gpoint_flavor_tmp);
    this->flavor        = Kokkos::create_mirror_view(DefaultDevice(), flavor_tmp);
    this->gpoint_flavor = Kokkos::create_mirror_view(DefaultDevice(), gpoint_flavor_tmp);
    Kokkos::deep_copy(this->flavor, flavor_tmp);
    Kokkos::deep_copy(this->gpoint_flavor, gpoint_flavor_tmp);

    // minimum, maximum reference temperature, pressure -- assumes low-to-high ordering
    //   for T, high-to-low ordering for p
    this->temp_ref_min  = temp_ref (0);
    this->temp_ref_max  = temp_ref (temp_ref.extent(0)-1);
    this->press_ref_min = press_ref(press_ref.extent(0)-1);
    this->press_ref_max = press_ref(0);

    // creates press_ref_log, temp_ref_delta
    this->press_ref_log_delta = (log(this->press_ref_min)-log(this->press_ref_max))/(press_ref.extent(0)-1);
    this->temp_ref_delta      = (this->temp_ref_max-this->temp_ref_min)/(temp_ref.extent(0)-1);

    // Which species are key in one or more bands?
    //   this->flavor is an index into this->gas_names
    this->is_key = bool1dk("is_key",this->get_ngas());
    // do j = 1, size(this%flavor, 2)
    //   do i = 1, size(this%flavor, 1) ! extents should be 2
    Kokkos::parallel_for( MDRangeP<2>( {0, 0}, {this->flavor.extent(0), this->flavor.extent(1)} ) , KOKKOS_LAMBDA (int i, int j) {
      if (this->flavor(i,j) != -1) { this->is_key(this->flavor(i,j)) = true; }
    });
  }

  // Initialize object based on data read from netCDF file however the user desires.
  //  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  // This interface is for the internal-sources object -- includes Plank functions and fractions
  void load(GasConcsK   const &available_gases,
            string1dv   const &gas_names,
            intHost3dk  const &key_species,
            intHost2dk  const &band2gpt,
            realHost2dk const &band_lims_wavenum,
            realHost1dk const &press_ref,
            real              press_ref_trop,
            realHost1dk const &temp_ref,
            real              temp_ref_p,
            real              temp_ref_t,
            realHost3dk const &vmr_ref,
            realHost4dk const &kmajor,
            realHost3dk const &kminor_lower,
            realHost3dk const &kminor_upper,
            string1dv   const &gas_minor,
            string1dv   const &identifier_minor,
            string1dv   const &minor_gases_lower,
            string1dv   const &minor_gases_upper,
            intHost2dk  const &minor_limits_gpt_lower,
            intHost2dk  const &minor_limits_gpt_upper,
            boolHost1dk const &minor_scales_with_density_lower,
            boolHost1dk const &minor_scales_with_density_upper,
            string1dv   const &scaling_gas_lower,
            string1dv   const &scaling_gas_upper,
            boolHost1dk const &scale_by_complement_lower,
            boolHost1dk const &scale_by_complement_upper,
            intHost1dk  const &kminor_start_lower,
            intHost1dk  const &kminor_start_upper,
            realHost2dk const &totplnk,
            realHost4dk const &planck_frac,
            realHost3dk const &rayl_lower,
            realHost3dk const &rayl_upper) {
    init_abs_coeffs(available_gases, gas_names, key_species, band2gpt, band_lims_wavenum, press_ref, temp_ref,
                    press_ref_trop, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
                    gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,
                    minor_limits_gpt_upper, minor_scales_with_density_lower, minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,
                    kminor_start_lower, kminor_start_upper, rayl_lower, rayl_upper);

    // Planck function tables
    this->totplnk = Kokkos::create_mirror_view(DefaultDevice(), totplnk);
    this->planck_frac = Kokkos::create_mirror_view(DefaultDevice(), planck_frac);
    Kokkos::deep_copy(this->totplnk, totplnk);
    Kokkos::deep_copy(this->planck_frac, planck_frac);
    // Temperature steps for Planck function interpolation
    //   Assumes that temperature minimum and max are the same for the absorption coefficient grid and the
    //   Planck grid and the Planck grid is equally spaced
    this->totplnk_delta = (this->temp_ref_max - this->temp_ref_min) / (this->totplnk.extent(0)-1);
  }



  // Initialize object based on data read from netCDF file however the user desires.
  //  Rayleigh scattering tables may or may not be present; this is indicated with allocation status
  // This interface is for the external-sources object -- includes TOA source function table
  void load(GasConcsK   const &available_gases,
            string1dv   const &gas_names,
            intHost3dk  const &key_species,
            intHost2dk  const &band2gpt,
            realHost2dk const &band_lims_wavenum,
            realHost1dk const &press_ref,
            real              press_ref_trop,
            realHost1dk const &temp_ref,
            real              temp_ref_p,
            real              temp_ref_t,
            realHost3dk const &vmr_ref,
            realHost4dk const &kmajor,
            realHost3dk const &kminor_lower,
            realHost3dk const &kminor_upper,
            string1dv   const &gas_minor,
            string1dv   const &identifier_minor,
            string1dv   const &minor_gases_lower,
            string1dv   const &minor_gases_upper,
            intHost2dk  const &minor_limits_gpt_lower,
            intHost2dk  const &minor_limits_gpt_upper,
            boolHost1dk const &minor_scales_with_density_lower,
            boolHost1dk const &minor_scales_with_density_upper,
            string1dv   const &scaling_gas_lower,
            string1dv   const &scaling_gas_upper,
            boolHost1dk const &scale_by_complement_lower,
            boolHost1dk const &scale_by_complement_upper,
            intHost1dk  const &kminor_start_lower,
            intHost1dk  const &kminor_start_upper,
            realHost1dk const &solar_src,
            realHost3dk const &rayl_lower,
            realHost3dk const &rayl_upper) {
    init_abs_coeffs(available_gases,  gas_names, key_species, band2gpt, band_lims_wavenum, press_ref, temp_ref,
                    press_ref_trop, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
                    gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,
                    minor_limits_gpt_upper, minor_scales_with_density_lower, minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,
                    kminor_start_lower, kminor_start_upper, rayl_lower, rayl_upper);

    // Solar source table init
    this->solar_src = Kokkos::create_mirror_view(DefaultDevice(), solar_src);
    Kokkos::deep_copy(this->solar_src, solar_src);
    this->totplnk_delta = 0.;
  }

  // Two functions to define array sizes needed by gas_optics()
  int get_ngas() const { return this->gas_names.size(); }

  // return the number of distinct major gas pairs in the spectral bands (referred to as
  // "flavors" - all bands have a flavor even if there is one or no major gas)
  int get_nflav() const { return this->flavor.extent(1); }

  string1dv get_gases() const { return this->gas_names; }

  // return the minimum pressure on the interpolation grids
  real get_press_min() const { return this->press_ref_min; }

  // return the maximum pressure on the interpolation grids
  real get_press_max() const { return this->press_ref_max; }

  // return the minimum temparature on the interpolation grids
  real get_temp_min()const  { return this->temp_ref_min; }

  // return the maximum temparature on the interpolation grids
  real get_temp_max() const { return this->temp_ref_max; }

  int get_neta() const { return this->kmajor.extent(1); }

  // return the number of pressures in reference profile
  //   absorption coefficient table is one bigger since a pressure is repeated in upper/lower atmos
  int get_npres() const { return this->kmajor.extent(2)-1; }

  int get_ntemp() const { return this->kmajor.extent(3); }

  // return the number of temperatures for Planck function
  int get_nPlanckTemp() const { return this->totplnk.extent(0); }


  // Function to define names of key and minor gases to be used by gas_optics().
  // The final list gases includes those that are defined in gas_optics_specification
  // and are provided in ty_gas_concs.
  string1dv get_minor_list(GasConcsK const &gas_desc, int ngas, string1dv const &name_spec) const {
    // List of minor gases to be used in gas_optics()
    string1dv rv;
    for (int igas=0 ; igas < this->get_ngas() ; igas++) {
      if (string_in_array(name_spec[igas], gas_desc.gas_name)) {
        rv.push_back(this->gas_names[igas]);
      }
    }
    return rv;
  }

  // return true if initialized for internal sources, false otherwise
  bool source_is_internal() const { return this->totplnk.is_allocated() && this->planck_frac.is_allocated(); }

  // return true if initialized for external sources, false otherwise
  bool source_is_external() const { return this->solar_src.is_allocated(); }

  // Ensure that every key gas required by the k-distribution is present in the gas concentration object
  void check_key_species_present(GasConcsK const &gas_desc) const {
    string1dv key_gas_names;
    auto is_key_h = Kokkos::create_mirror_view(this->is_key);
    Kokkos::deep_copy(is_key_h, this->is_key);
    for (auto i = 0; i < is_key_h.extent(0); ++i) {
      if (is_key_h(i)) {
        key_gas_names.push_back(this->gas_names[i]);
      }
    }
    for (auto igas=0 ; igas < key_gas_names.size() ; igas++) {
      if (! string_in_array(key_gas_names[igas], gas_desc.gas_name)) {
        stoprun("gas required by k-distribution is not present in the GasConcs object");
      }
    }
  }

  // Compute gas optical depth and Planck source functions, given temperature, pressure, and composition
  template <class T>
  void gas_optics(const int ncol, const int nlay,
                  bool top_at_1, real2dk const &play, real2dk const &plev, real2dk const &tlay, real1dk const &tsfc,
                  GasConcsK const &gas_desc, T &optical_props, SourceFuncLWK &sources,
                  real2dk const &col_dry=real2dk(), real2dk const &tlev=real2dk()) {
    int ngpt  = this->get_ngpt();
    int nband = this->get_nband();
    // Interpolation coefficients for use in source function
    int2dk  jtemp ("jtemp"                         ,ncol,nlay);
    int2dk  jpress("jpress"                        ,ncol,nlay);
    bool2dk tropo ("tropo"                         ,ncol,nlay);
    real6dk fmajor("fmajor",2,2,2,this->get_nflav(),ncol,nlay);
    int4dk  jeta  ("jeta"  ,2    ,this->get_nflav(),ncol,nlay);
    // Gas optics
    compute_gas_taus(top_at_1, ncol, nlay, ngpt, nband, play, plev, tlay, gas_desc, optical_props, jtemp, jpress,
                     jeta, tropo, fmajor, col_dry);

    // External source -- check arrays sizes and values
    // input data sizes and values
    if (tsfc.extent(0) != ncol) { stoprun("gas_optics(): array tsfc has wrong size"); }
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if (any(tsfc < this->temp_ref_min) || any(tsfc > this->temp_ref_max)) {
        stoprun("gas_optics(): array tsfc has values outside range");
      }
    #endif

    if (tlev.is_allocated()) {
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (any(tlev < this->temp_ref_min) || any(tlev > this->temp_ref_max)) {
          stoprun("gas_optics(): array tlev has values outside range");
        }
      #endif
    }

    // output extents
    if (sources.get_ncol() != ncol || sources.get_nlay() != nlay || sources.get_ngpt() != ngpt) {
      stoprun("gas_optics%gas_optics: source function arrays inconsistently sized");
    }

    // Interpolate source function
    this->source(top_at_1, ncol, nlay, nband, ngpt, play, plev, tlay, tsfc, jtemp, jpress, jeta, tropo, fmajor, sources, tlev);
  }

  // Compute gas optical depth given temperature, pressure, and composition
  template <class T>
  void gas_optics(const int ncol, const int nlay,
                  bool top_at_1, real2dk const &play, real2dk const &plev, real2dk const &tlay, GasConcsK const &gas_desc,
                  T &optical_props, real2dk &toa_src, real2dk const &col_dry=real2dk()) {
    int ngpt  = this->get_ngpt();
    int nband = this->get_nband();
    int ngas  = this->get_ngas();
    int nflav = get_nflav();

    // Interpolation coefficients for use in source function
    int2dk  jtemp ("jtemp"                         ,ncol,nlay);
    int2dk  jpress("jpress"                        ,ncol,nlay);
    bool2dk tropo ("tropo"                         ,ncol,nlay);
    real6dk fmajor("fmajor",2,2,2,this->get_nflav(),ncol,nlay);
    int4dk  jeta  ("jeta  ",2    ,this->get_nflav(),ncol,nlay);
    // Gas optics
    compute_gas_taus(top_at_1, ncol, nlay, ngpt, nband, play, plev, tlay, gas_desc, optical_props, jtemp, jpress, jeta,
                     tropo, fmajor, col_dry);

    // External source function is constant
    if (toa_src.extent(0) != ncol || toa_src.extent(1) != ngpt) { stoprun("gas_optics(): array toa_src has wrong size"); }

    Kokkos::parallel_for( MDRangeP<2>({0,0}, {ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      toa_src(icol,igpt) = this->solar_src(igpt);
    });
  }

  // Returns optical properties and interpolation coefficients
  template <class T>
  void compute_gas_taus(bool top_at_1, int ncol, int nlay, int ngpt, int nband, real2dk const &play, real2dk const &plev, real2dk const &tlay,
                        GasConcsK const &gas_desc, T &optical_props, int2dk const &jtemp, int2dk const &jpress, int4dk const &jeta,
                        bool2dk const &tropo, real6dk const &fmajor, real2dk const &col_dry=real2dk() ) {
    // Number of molecules per cm^2
    real3dk tau         ("tau"         ,ngpt,nlay,ncol);
    real3dk tau_rayleigh("tau_rayleigh",ngpt,nlay,ncol);
    // Interpolation variables used in major gas but not elsewhere, so don't need exporting
    real3dk vmr         ("vmr"         ,ncol,nlay,this->get_ngas());
    realOff3dk col_gas  ("col_gas"     ,std::make_pair(0, ncol-1), std::make_pair(0, nlay-1), std::make_pair(-1, this->get_ngas()-1));
    real4dk col_mix     ("col_mix"     ,2,this->get_nflav(),ncol,nlay); // combination of major species's column amounts
                                                                       // index(1) : reference temperature level
                                                                       // index(2) : flavor
                                                                       // index(3) : layer
    real5dk fminor      ("fminor"      ,2,2,this->get_nflav(),ncol,nlay); // interpolation fractions for minor species
                                                                         // index(1) : reference eta level (temperature dependent)
                                                                         // index(2) : reference temperature level
                                                                         // index(3) : flavor
                                                                         // index(4) : layer
    // Error checking
    // Check for initialization
    if (! this->is_initialized()) { stoprun("ERROR: spectral configuration not loaded"); }
    // Check for presence of key species in ty_gas_concs; return error if any key species are not present
    this->check_key_species_present(gas_desc);
    #ifdef RRTMGP_EXPENSIVE_CHECKS
      if ( any(play < this->press_ref_min) || any(play > this->press_ref_max) ) {
        stoprun("gas_optics(): array play has values outside range");
      }
      if ( any(plev < this->press_ref_min) || any(plev > this->press_ref_max) ) {
        stoprun("gas_optics(): array plev has values outside range");
      }
      if ( any(tlay < this->temp_ref_min) || any(tlay > this->temp_ref_max) ) {
        stoprun("gas_optics(): array tlay has values outside range");
      }
    #endif
    if (col_dry.is_allocated()) {
      if (col_dry.extent(0) != ncol || col_dry.extent(1) != nlay) { stoprun("gas_optics(): array col_dry has wrong size"); }
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (any(col_dry < 0.)) { stoprun("gas_optics(): array col_dry has values outside range"); }
      #endif
    }

    bool use_rayl = this->krayl.is_allocated();

    int ngas  = this->get_ngas();
    int nflav = this->get_nflav();
    int neta  = this->get_neta();
    int npres = this->get_npres();
    int ntemp = this->get_ntemp();
    // number of minor contributors, total num absorption coeffs
    int nminorlower  = this->minor_scales_with_density_lower.extent(0);
    int nminorklower = this->kminor_lower.extent(0);
    int nminorupper  = this->minor_scales_with_density_upper.extent(0);
    int nminorkupper = this->kminor_upper.extent(0);
    // Fill out the array of volume mixing ratios
    for (int igas = 0 ; igas < ngas ; igas++) {
      // Get vmr if  gas is provided in ty_gas_concs
      for (size_t igas2 = 0 ; igas2 < gas_desc.gas_name.size() ; igas2++) {
        if ( lower_case(this->gas_names[igas]) == lower_case(gas_desc.gas_name[igas2]) ) {
          real2dk vmr_slice = Kokkos::subview(vmr, Kokkos::ALL, Kokkos::ALL, igas);
          gas_desc.get_vmr(this->gas_names[igas], vmr_slice);
        }
      }
    }

    // Compute dry air column amounts (number of molecule per cm^2) if user hasn't provided them
    int idx_h2o = string_loc_in_array("h2o", this->gas_names);
    real2dk col_dry_wk;
    if (col_dry.is_allocated()) {
      col_dry_wk = col_dry;
    } else {
      real2dk col_dry_arr = this->get_col_dry(Kokkos::subview(vmr, Kokkos::ALL, Kokkos::ALL, idx_h2o),plev); // dry air column amounts computation
      col_dry_wk = col_dry_arr;
    }
    // compute column gas amounts [molec/cm^2]
    // do ilay = 1, nlay
    //   do icol = 1, ncol
    Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay,ncol}) , KOKKOS_LAMBDA (int ilay, int icol) {
      col_gas(icol,ilay,-1) = col_dry_wk(icol,ilay);
    });
    // do igas = 1, ngas
    //   do ilay = 1, nlay
    //     do icol = 1, ncol
    Kokkos::parallel_for( MDRangeP<3>({0,0,0}, {ngas,nlay,ncol}) , KOKKOS_LAMBDA (int igas, int ilay, int icol) {
      col_gas(icol,ilay,igas) = vmr(icol,ilay,igas) * col_dry_wk(icol,ilay);
    });
    // ---- calculate gas optical depths ----
    Kokkos::deep_copy(tau, 0);

    interpolation(ncol, nlay, ngas, nflav, neta, npres, ntemp, this->flavor, this->press_ref_log, this->temp_ref,
                  this->press_ref_log_delta, this->temp_ref_min, this->temp_ref_delta, this->press_ref_trop_log,
                  this->vmr_ref, play, tlay, col_gas, jtemp, fmajor, fminor, col_mix, tropo, jeta, jpress);

    compute_tau_absorption(this->max_gpt_diff_lower, this->max_gpt_diff_upper, ncol, nlay, nband, ngpt, ngas, nflav, neta, npres, ntemp, nminorlower, nminorklower,
                           nminorupper, nminorkupper, idx_h2o, this->gpoint_flavor, this->get_band_lims_gpoint(),
                           this->kmajor, this->kminor_lower, this->kminor_upper, this->minor_limits_gpt_lower,
                           this->minor_limits_gpt_upper, this->minor_scales_with_density_lower,
                           this->minor_scales_with_density_upper, this->scale_by_complement_lower,
                           this->scale_by_complement_upper, this->idx_minor_lower, this->idx_minor_upper,
                           this->idx_minor_scaling_lower, this->idx_minor_scaling_upper, this->kminor_start_lower,
                           this->kminor_start_upper, tropo, col_mix, fmajor, fminor, play, tlay, col_gas,
                           jeta, jtemp, jpress, tau, top_at_1);

    if (this->krayl.is_allocated()) {
      compute_tau_rayleigh( ncol, nlay, nband, ngpt, ngas, nflav, neta, npres, ntemp, this->gpoint_flavor,
                            this->get_band_lims_gpoint(), this->krayl, idx_h2o, col_dry_wk, col_gas,
                            fminor, jeta, tropo, jtemp, tau_rayleigh);
    }
    combine_and_reorder(tau, tau_rayleigh, this->krayl.is_allocated(), optical_props);
  }

  // Compute Planck source functions at layer centers and levels
  void source(bool top_at_1, int ncol, int nlay, int nbnd, int ngpt, real2dk const &play, real2dk const &plev, real2dk const &tlay,
              real1dk const &tsfc, int2dk const &jtemp, int2dk const &jpress, int4dk const &jeta, bool2dk const &tropo,
              real6dk const &fmajor, SourceFuncLWK &sources, real2dk const &tlev=real2dk()) {
    real3dk lay_source_t    ("lay_source_t    ",ngpt,nlay,ncol);
    real3dk lev_source_inc_t("lev_source_inc_t",ngpt,nlay,ncol);
    real3dk lev_source_dec_t("lev_source_dec_t",ngpt,nlay,ncol);
    real2dk sfc_source_t    ("sfc_source_t    ",ngpt     ,ncol);
    // Variables for temperature at layer edges [K] (ncol, nlay+1)
    real2dk tlev_arr("tlev_arr",ncol,nlay+1);

    // Source function needs temperature at interfaces/levels and at layer centers
    real2dk tlev_wk;
    if (tlev.is_allocated()) {
      //   Users might have provided these
      tlev_wk = tlev;
    } else {
      tlev_wk = real2dk("tlev_wk",ncol,nlay+1);
      // Interpolate temperature to levels if not provided
      //   Interpolation and extrapolation at boundaries is weighted by pressure
      // do ilay = 1, nlay+1
      //   do icol = 1, ncol
      Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlay+1,ncol}) , KOKKOS_LAMBDA (int ilay, int icol) {
        if (ilay == 0) {
          tlev_wk(icol,0) = tlay(icol,0) + (plev(icol,0)-play(icol,0))*(tlay(icol,1)-tlay(icol,0)) / (play(icol,1)-play(icol,0));
        }
        else if (ilay == nlay) {
          tlev_wk(icol,ilay) = tlay(icol,ilay-1) + (plev(icol,ilay)-play(icol,ilay-1))*(tlay(icol,ilay-1)-tlay(icol,ilay)) /
            (play(icol,nlay)-play(icol,ilay));
        }
        else {
          tlev_wk(icol,ilay) = ( play(icol,ilay-1)*tlay(icol,ilay-1)*(plev(icol,ilay  )-play(icol,ilay))  +
                                 play(icol,ilay  )*tlay(icol,ilay  )*(play(icol,ilay-1)-plev(icol,ilay)) ) /
            (plev(icol,ilay)*(play(icol,ilay-1) - play(icol,ilay)));
        }
      });
    }
    // Compute internal (Planck) source functions at layers and levels,
    //  which depend on mapping from spectral space that creates k-distribution.
    int nlayTmp = conv::merge( nlay-1 , 0 , top_at_1 );
    compute_Planck_source(ncol, nlay, nbnd, ngpt, this->get_nflav(), this->get_neta(), this->get_npres(), this->get_ntemp(),
                          this->get_nPlanckTemp(), tlay, tlev_wk, tsfc, nlayTmp, fmajor, jeta, tropo, jtemp, jpress,
                          this->get_gpoint_bands(), this->get_band_lims_gpoint(), this->planck_frac, this->temp_ref_min,
                          this->totplnk_delta, this->totplnk, this->gpoint_flavor, sfc_source_t, lay_source_t, lev_source_inc_t,
                          lev_source_dec_t);

    auto &sources_sfc_source = sources.sfc_source;
    // do igpt = 1, ngpt
    //   do icol = 1, ncol
    Kokkos::parallel_for( MDRangeP<2>({0,0}, {ngpt,ncol}) , KOKKOS_LAMBDA (int igpt, int icol) {
      sources_sfc_source(icol,igpt) = sfc_source_t(igpt,icol);
    });
    reorder123x321(ngpt, nlay, ncol, lay_source_t    , sources.lay_source    );
    reorder123x321(ngpt, nlay, ncol, lev_source_inc_t, sources.lev_source_inc);
    reorder123x321(ngpt, nlay, ncol, lev_source_dec_t, sources.lev_source_dec);
  }

  // Utility function, provided for user convenience
  // computes column amounts of dry air using hydrostatic equation
  real2dk get_col_dry(real2dk const &vmr_h2o, real2dk const &plev, real1dk const &latitude=real1dk()) {
    // first and second term of Helmert formula
    real constexpr helmert1 = 9.80665;
    real constexpr helmert2 = 0.02586;
    int ncol = plev.extent(0);
    int nlev = plev.extent(1);
    real1dk g0("g0",plev.extent(0));
    if (latitude.is_allocated()) {
      // A purely OpenACC implementation would probably compute g0 within the kernel below
      // do icol = 1, ncol
      Kokkos::parallel_for( ncol , KOKKOS_LAMBDA (int icol) {
        g0(icol) = helmert1 - helmert2 * cos(2.0 * M_PI * latitude(icol) / 180.0); // acceleration due to gravity [m/s^2]
      });
    } else {
      // do icol = 1, ncol
      const auto grav = ::grav;
      Kokkos::parallel_for( ncol, KOKKOS_LAMBDA (int icol) {
        g0(icol) = grav;
      });
    }

    real2dk col_dry("col_dry",plev.extent(0),plev.extent(1)-1);
    // do ilev = 1, nlev-1
    //   do icol = 1, ncol
    const auto m_dry = ::m_dry;
    Kokkos::parallel_for( MDRangeP<2>({0,0}, {nlev-1,ncol}) , KOKKOS_LAMBDA (int ilev , int icol) {
      real delta_plev = abs(plev(icol,ilev) - plev(icol,ilev+1));
      // Get average mass of moist air per mole of moist air
      real fact = 1. / (1.+vmr_h2o(icol,ilev));
      real m_air = (m_dry + m_h2o * vmr_h2o(icol,ilev)) * fact;
      col_dry(icol,ilev) = 10. * delta_plev * avogad * fact/(1000.*m_air*100.*g0(icol));
    });
    return col_dry;
  }

 // Utility function to combine optical depths from gas absorption and Rayleigh scattering
 //   (and reorder them for convenience, while we're at it)
 void combine_and_reorder(real3dk const &tau, real3dk const &tau_rayleigh, bool has_rayleigh, OpticalProps1sclK &optical_props) {
    int ncol = tau.extent(2);
    int nlay = tau.extent(1);
    int ngpt = tau.extent(0);
    reorder123x321(ngpt, nlay, ncol, tau, optical_props.tau);
  }

 // Utility function to combine optical depths from gas absorption and Rayleigh scattering
 //   (and reorder them for convenience, while we're at it)
 void combine_and_reorder(real3dk const &tau, real3dk const &tau_rayleigh, bool has_rayleigh, OpticalProps2strK &optical_props) {
    int ncol = tau.extent(2);
    int nlay = tau.extent(1);
    int ngpt = tau.extent(0);
    if (has_rayleigh) {
      // combine optical depth and rayleigh scattering
      combine_and_reorder_2str(ncol, nlay, ngpt, tau, tau_rayleigh, optical_props.tau, optical_props.ssa, optical_props.g);
    } else {
      // index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
      reorder123x321(ngpt, nlay, ncol, tau, optical_props.tau);
      Kokkos::deep_copy(optical_props.ssa, 0);
      Kokkos::deep_copy(optical_props.g,   0);
    }
  }

  void print_norms() const {
                                                      std::cout << "name                                  : " << std::setw(20) << name                                   << "\n";
                                                      std::cout << "totplnk_delta                         : " << std::setw(20) << totplnk_delta                          << "\n";
                                                      std::cout << "press_ref_min                         : " << std::setw(20) << press_ref_min                          << "\n";
                                                      std::cout << "press_ref_max                         : " << std::setw(20) << press_ref_max                          << "\n";
                                                      std::cout << "temp_ref_min                          : " << std::setw(20) << temp_ref_min                           << "\n";
                                                      std::cout << "temp_ref_max                          : " << std::setw(20) << temp_ref_max                           << "\n";
                                                      std::cout << "press_ref_log_delta                   : " << std::setw(20) << press_ref_log_delta                    << "\n";
                                                      std::cout << "temp_ref_delta                        : " << std::setw(20) << temp_ref_delta                         << "\n";
                                                      std::cout << "press_ref_trop_log                    : " << std::setw(20) << press_ref_trop_log                     << "\n";
    //if (gas_names.is_allocated()                      ) { std::cout << "gas_names                             : " << std::setw(20) << gas_names                              << "\n"; }
    if (band2gpt.is_allocated()                       ) { std::cout << "sum(band2gpt     )                    : " << std::setw(20) << conv::sum(band2gpt     )                     << "\n"; }
    if (gpt2band.is_allocated()                       ) { std::cout << "sum(gpt2band     )                    : " << std::setw(20) << conv::sum(gpt2band     )                     << "\n"; }
    if (band_lims_wvn.is_allocated()                  ) { std::cout << "sum(band_lims_wvn)                    : " << std::setw(20) << conv::sum(band_lims_wvn)                     << "\n"; }
    if (press_ref.is_allocated()                      ) { std::cout << "sum(press_ref    )                    : " << std::setw(20) << conv::sum(press_ref    )                     << "\n"; }
    if (press_ref_log.is_allocated()                  ) { std::cout << "sum(press_ref_log)                    : " << std::setw(20) << conv::sum(press_ref_log)                     << "\n"; }
    if (temp_ref.is_allocated()                       ) { std::cout << "sum(temp_ref     )                    : " << std::setw(20) << conv::sum(temp_ref     )                     << "\n"; }
    if (vmr_ref.is_allocated()                        ) { std::cout << "sum(vmr_ref                )          : " << std::setw(20) << conv::sum(vmr_ref                )           << "\n"; }
    if (flavor.is_allocated()                         ) { std::cout << "sum(flavor                 )          : " << std::setw(20) << conv::sum(flavor                 )           << "\n"; }
    if (gpoint_flavor.is_allocated()                  ) { std::cout << "sum(gpoint_flavor          )          : " << std::setw(20) << conv::sum(gpoint_flavor          )           << "\n"; }
    if (kmajor.is_allocated()                         ) { std::cout << "sum(kmajor                 )          : " << std::setw(20) << conv::sum(kmajor                 )           << "\n"; }
    if (minor_limits_gpt_lower.is_allocated()         ) { std::cout << "sum(minor_limits_gpt_lower )          : " << std::setw(20) << conv::sum(minor_limits_gpt_lower )           << "\n"; }
    if (minor_limits_gpt_upper.is_allocated()         ) { std::cout << "sum(minor_limits_gpt_upper )          : " << std::setw(20) << conv::sum(minor_limits_gpt_upper )           << "\n"; }
    if (idx_minor_lower.is_allocated()                ) { std::cout << "sum(idx_minor_lower        )          : " << std::setw(20) << conv::sum(idx_minor_lower        )           << "\n"; }
    if (idx_minor_upper.is_allocated()                ) { std::cout << "sum(idx_minor_upper        )          : " << std::setw(20) << conv::sum(idx_minor_upper        )           << "\n"; }
    if (idx_minor_scaling_lower.is_allocated()        ) { std::cout << "sum(idx_minor_scaling_lower)          : " << std::setw(20) << conv::sum(idx_minor_scaling_lower)           << "\n"; }
    if (idx_minor_scaling_upper.is_allocated()        ) { std::cout << "sum(idx_minor_scaling_upper)          : " << std::setw(20) << conv::sum(idx_minor_scaling_upper)           << "\n"; }
    if (kminor_start_lower.is_allocated()             ) { std::cout << "sum(kminor_start_lower     )          : " << std::setw(20) << conv::sum(kminor_start_lower     )           << "\n"; }
    if (kminor_start_upper.is_allocated()             ) { std::cout << "sum(kminor_start_upper     )          : " << std::setw(20) << conv::sum(kminor_start_upper     )           << "\n"; }
    if (kminor_lower.is_allocated()                   ) { std::cout << "sum(kminor_lower           )          : " << std::setw(20) << conv::sum(kminor_lower           )           << "\n"; }
    if (kminor_upper.is_allocated()                   ) { std::cout << "sum(kminor_upper           )          : " << std::setw(20) << conv::sum(kminor_upper           )           << "\n"; }
    if (krayl.is_allocated()                          ) { std::cout << "sum(krayl                  )          : " << std::setw(20) << conv::sum(krayl                  )           << "\n"; }
    if (planck_frac.is_allocated()                    ) { std::cout << "sum(planck_frac            )          : " << std::setw(20) << conv::sum(planck_frac            )           << "\n"; }
    if (totplnk.is_allocated()                        ) { std::cout << "sum(totplnk                )          : " << std::setw(20) << conv::sum(totplnk                )           << "\n"; }
    if (solar_src.is_allocated()                      ) { std::cout << "sum(solar_src              )          : " << std::setw(20) << conv::sum(solar_src              )           << "\n"; }
    if (minor_scales_with_density_lower.is_allocated()) { std::cout << "count(minor_scales_with_density_lower): " << std::setw(20) << conv::sum(minor_scales_with_density_lower) << "\n"; }
    if (minor_scales_with_density_upper.is_allocated()) { std::cout << "count(minor_scales_with_density_upper): " << std::setw(20) << conv::sum(minor_scales_with_density_upper) << "\n"; }
    if (scale_by_complement_lower.is_allocated()      ) { std::cout << "count(scale_by_complement_lower      ): " << std::setw(20) << conv::sum(scale_by_complement_lower      ) << "\n"; }
    if (scale_by_complement_upper.is_allocated()      ) { std::cout << "count(scale_by_complement_upper      ): " << std::setw(20) << conv::sum(scale_by_complement_upper      ) << "\n"; }
    if (is_key.is_allocated()                         ) { std::cout << "count(is_key                         ): " << std::setw(20) << conv::sum(is_key                         ) << "\n"; }
  }

#ifdef RRTMGP_ENABLE_YAKL
  void validate_kokkos(const GasOpticsRRTMGP& orig) const
  {
    OpticalPropsK::validate_kokkos(orig);

    conv::compare_yakl_to_kokkos(orig.press_ref, press_ref);
    conv::compare_yakl_to_kokkos(orig.press_ref_log, press_ref_log);
    conv::compare_yakl_to_kokkos(orig.temp_ref, temp_ref);

    RRT_REQUIRE(orig.press_ref_min == press_ref_min, "Bad press_ref_min");
    RRT_REQUIRE(orig.press_ref_max == press_ref_max, "Bad press_ref_max");
    RRT_REQUIRE(orig.temp_ref_min == temp_ref_min, "Bad temp_ref_min");
    RRT_REQUIRE(orig.temp_ref_max == temp_ref_max, "Bad temp_ref_max");
    RRT_REQUIRE(orig.press_ref_log_delta == press_ref_log_delta, "Bad press_ref_log_delta");
    RRT_REQUIRE(orig.temp_ref_delta == temp_ref_delta, "Bad temp_ref_delta");
    RRT_REQUIRE(orig.press_ref_trop_log == press_ref_trop_log, "Bad press_ref_trop_log");
    RRT_REQUIRE(orig.totplnk_delta == totplnk_delta, "Bad totplnk_delta");

    RRT_REQUIRE(orig.max_gpt_diff_lower == max_gpt_diff_lower, "Bad max_gpt_diff_lower");
    RRT_REQUIRE(orig.max_gpt_diff_upper == max_gpt_diff_upper, "Bad max_gpt_diff_upper");

    conv::compare_yakl_to_kokkos_str(orig.gas_names, gas_names);
    conv::compare_yakl_to_kokkos(orig.vmr_ref, vmr_ref);

    conv::compare_yakl_to_kokkos(orig.flavor, flavor, true);
    conv::compare_yakl_to_kokkos(orig.gpoint_flavor, gpoint_flavor, true /*idx data*/);

    conv::compare_yakl_to_kokkos(orig.kmajor, kmajor);

    conv::compare_yakl_to_kokkos(orig.minor_limits_gpt_lower, minor_limits_gpt_lower, true);
    conv::compare_yakl_to_kokkos(orig.minor_limits_gpt_upper, minor_limits_gpt_upper, true);

    conv::compare_yakl_to_kokkos(orig.minor_scales_with_density_lower, minor_scales_with_density_lower);
    conv::compare_yakl_to_kokkos(orig.minor_scales_with_density_upper, minor_scales_with_density_upper);
    conv::compare_yakl_to_kokkos(orig.scale_by_complement_lower, scale_by_complement_lower);
    conv::compare_yakl_to_kokkos(orig.scale_by_complement_upper, scale_by_complement_upper);
    conv::compare_yakl_to_kokkos(orig.idx_minor_lower, idx_minor_lower, true);
    conv::compare_yakl_to_kokkos(orig.idx_minor_upper, idx_minor_upper, true);
    conv::compare_yakl_to_kokkos(orig.idx_minor_scaling_lower, idx_minor_scaling_lower, true);
    conv::compare_yakl_to_kokkos(orig.idx_minor_scaling_upper, idx_minor_scaling_upper, true);

    conv::compare_yakl_to_kokkos(orig.kminor_start_lower, kminor_start_lower, true);
    conv::compare_yakl_to_kokkos(orig.kminor_start_upper, kminor_start_upper, true);

    conv::compare_yakl_to_kokkos(orig.kminor_lower, kminor_lower);
    conv::compare_yakl_to_kokkos(orig.kminor_upper, kminor_upper);

    conv::compare_yakl_to_kokkos(orig.krayl, krayl);
    conv::compare_yakl_to_kokkos(orig.planck_frac, planck_frac);

    conv::compare_yakl_to_kokkos(orig.totplnk, totplnk);

    conv::compare_yakl_to_kokkos(orig.solar_src, solar_src);

    conv::compare_yakl_to_kokkos(orig.is_key, is_key);
  }
#endif

};
#endif
