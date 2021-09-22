
#pragma once

#include "mo_optical_props.h"
#include "mo_source_functions.h"
#include "mo_rrtmgp_util_string.h"
#include "mo_gas_optics_kernels.h"
#include "mo_rrtmgp_constants.h"
#include "mo_rrtmgp_util_reorder.h"
#include "mo_gas_concentrations.h"

using yakl::intrinsics::count;
using yakl::intrinsics::pack;
using yakl::intrinsics::lbound;
using yakl::intrinsics::ubound;
using yakl::COLON;

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
    idx_minor_scaling_atm = intHost1d("idx_minor_scaling_atm",size(scaling_gas_atm,1));
    for (int imnr=1 ; imnr <= size(scaling_gas_atm,1) ; imnr++) {
      // This will be -1 if there's no interacting gas
      idx_minor_scaling_atm(imnr) = string_loc_in_array(scaling_gas_atm(imnr), gas_names);
    }
  }



  void create_key_species_reduce(string1d const &gas_names, string1d const &gas_names_red, intHost3d const &key_species,
                                 intHost3d &key_species_red, boolHost1d &key_species_present_init) {
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
            key_species_red(ip,ia,it) = key_species(ip,ia,it);
          }
        }
      }
    }
  }



  // Create flavor list
  // An unordered array of extent (2,:) containing all possible pairs of key species used in either upper or lower atmos
  void create_flavor(intHost3d const &key_species, intHost2d &flavor) {
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
    int ngpt = size(gpt2band,1);
    gpoint_flavor = intHost2d("gpoint_flavor",2,ngpt);
    for (int igpt=1 ; igpt <= ngpt ; igpt++) {
      for (int iatm = 1 ; iatm <= 2 ; iatm++) {
        int key_species_pair2flavor = -1;
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
      this->max_gpt_diff_lower = max( this->max_gpt_diff_lower , minor_limits_gpt_lower_red(2,i) - minor_limits_gpt_lower_red(1,i) );
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
      this->max_gpt_diff_upper = max( this->max_gpt_diff_upper , minor_limits_gpt_upper_red(2,i) - minor_limits_gpt_upper_red(1,i) );
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
    auto &press_ref_loc     = this->press_ref;
    auto &press_ref_log_loc = this->press_ref_log;
    // Running a kernel because it's more convenient in this case
    parallel_for( Bounds<1>( size(this->press_ref,1) ) , YAKL_LAMBDA (int i) {
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
    auto &is_key_loc = this->is_key;
    auto &flavor_loc = this->flavor;
    parallel_for( Bounds<1>( this->get_ngas() ) , YAKL_LAMBDA (int i) {
      is_key_loc(i) = false;
    });
    // do j = 1, size(this%flavor, 2)
    //   do i = 1, size(this%flavor, 1) ! extents should be 2
    parallel_for( Bounds<2>( size(this->flavor,2) , size(this->flavor,1) ) , YAKL_LAMBDA (int j, int i) {
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

    init_abs_coeffs(available_gases,  gas_names, key_species, band2gpt, band_lims_wavenum, press_ref, temp_ref,
                    press_ref_trop, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper,
                    gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, minor_limits_gpt_lower,
                    minor_limits_gpt_upper, minor_scales_with_density_lower, minor_scales_with_density_upper,
                    scaling_gas_lower, scaling_gas_upper, scale_by_complement_lower, scale_by_complement_upper,
                    kminor_start_lower, kminor_start_upper, rayl_lower, rayl_upper);
    
    // Solar source table init
    this->solar_src = real1d("solar_src",size(solar_src,1));
    solar_src.deep_copy_to(this->solar_src);
  }



  // Two functions to define array sizes needed by gas_optics()
  int get_ngas() const { return size(this->gas_names,1); }



  // return the number of distinct major gas pairs in the spectral bands (referred to as
  // "flavors" - all bands have a flavor even if there is one or no major gas)
  int get_nflav() const { return size(this->flavor,2); }



  string1d get_gases() const { return this->gas_names; }



  // return the minimum pressure on the interpolation grids
  real get_press_min() const { return this->press_ref_min; }



  // return the maximum pressure on the interpolation grids
  real get_press_max() const { return this->press_ref_max; }



  // return the minimum temparature on the interpolation grids
  real get_temp_min()const  { return this->temp_ref_min; }



  // return the maximum temparature on the interpolation grids
  real get_temp_max() const { return this->temp_ref_max; }



  int get_neta() const { return size(this->kmajor,2); }



  // return the number of pressures in reference profile
  //   absorption coefficient table is one bigger since a pressure is repeated in upper/lower atmos
  int get_npres() const { return size(this->kmajor,3)-1; }



  int get_ntemp() const { return size(this->kmajor,4); }



  // return the number of temperatures for Planck function
  int get_nPlanckTemp() const { return size(this->totplnk,1); }



  // Function to define names of key and minor gases to be used by gas_optics().
  // The final list gases includes those that are defined in gas_optics_specification
  // and are provided in ty_gas_concs.
  string1d get_minor_list(GasConcs const &gas_desc, int ngas, string1d const &name_spec) const {
    // List of minor gases to be used in gas_optics()
    boolHost1d gas_is_present("gas_is_present",size(name_spec,1));
    for (int igas=1 ; igas <= this->get_ngas() ; igas++) {
      gas_is_present(igas) = string_in_array(name_spec(igas), gas_desc.gas_name);
    }
    return pack(this->gas_names, gas_is_present);
  }



  // return true if initialized for internal sources, false otherwise
  bool source_is_internal() const { return allocated(this->totplnk) && allocated(this->planck_frac); }



  // return true if initialized for external sources, false otherwise
  bool source_is_external() const { return allocated(this->solar_src); }



  // Ensure that every key gas required by the k-distribution is present in the gas concentration object
  void check_key_species_present(GasConcs const &gas_desc) const {
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
      if (anyLT(tsfc,this->temp_ref_min) || anyGT(tsfc,this->temp_ref_max)) {
        stoprun("gas_optics(): array tsfc has values outside range");
      }
    #endif

    if (allocated(tlev)) {
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (anyLT(tlev,this->temp_ref_min) || anyGT(tlev,this->temp_ref_max)) {
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

    auto &solar_src_loc = this->solar_src;
    parallel_for( Bounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      toa_src(icol,igpt) = solar_src_loc(igpt);
    });
  }



  // Returns optical properties and interpolation coefficients
  template <class T>
  void compute_gas_taus(bool top_at_1, int ncol, int nlay, int ngpt, int nband, real2d const &play, real2d const &plev, real2d const &tlay,
                        GasConcs const &gas_desc, T &optical_props, int2d &jtemp, int2d &jpress, int4d &jeta,
                        bool2d &tropo, real6d &fmajor, real2d const &col_dry=real2d() ) {
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
      if ( anyLT(play,this->press_ref_min) || anyGT(play,this->press_ref_max) ) {
        stoprun("gas_optics(): array play has values outside range");
      }
      if ( anyLT(plev,this->press_ref_min) || anyGT(plev,this->press_ref_max) ) {
        stoprun("gas_optics(): array plev has values outside range");
      }
      if ( anyLT(tlay,this->temp_ref_min) || anyGT(tlay,this->temp_ref_max) ) {
        stoprun("gas_optics(): array tlay has values outside range");
      }
    #endif
    if (allocated(col_dry)) {
      if (size(col_dry,1) != ncol || size(col_dry,2) != nlay) { stoprun("gas_optics(): array col_dry has wrong size"); }
      #ifdef RRTMGP_EXPENSIVE_CHECKS
        if (anyLT(col_dry,0._wp)) { stoprun("gas_optics(): array col_dry has values outside range"); }
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
      for (int igas2 = lbound(gas_desc.gas_name,1) ; igas2 <= ubound(gas_desc.gas_name,1) ; igas2++) {
        if ( lower_case(this->gas_names(igas)) == lower_case(gas_desc.gas_name(igas2)) ) {
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
    parallel_for( Bounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
      col_gas(icol,ilay,0) = col_dry_wk(icol,ilay);
    });
    // do igas = 1, ngas
    //   do ilay = 1, nlay
    //     do icol = 1, ncol
    parallel_for( Bounds<3>(ngas,nlay,ncol) , YAKL_LAMBDA (int igas, int ilay, int icol) {
      col_gas(icol,ilay,igas) = vmr(icol,ilay,igas) * col_dry_wk(icol,ilay);
    });
    // ---- calculate gas optical depths ----
    memset(tau , 0._wp);

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
      parallel_for( Bounds<2>(nlay+1,ncol) , YAKL_LAMBDA (int ilay, int icol) {
        if (ilay == 1) {
          tlev_wk(icol,1) = tlay(icol,1) + (plev(icol,1)-play(icol,1))*(tlay(icol,2)-tlay(icol,1)) / (play(icol,2)-play(icol,1));
        }
        tlev_wk(icol,ilay) = ( play(icol,ilay-1)*tlay(icol,ilay-1)*(plev(icol,ilay  )-play(icol,ilay))  +
                               play(icol,ilay  )*tlay(icol,ilay  )*(play(icol,ilay-1)-plev(icol,ilay)) ) /
                             (plev(icol,ilay)*(play(icol,ilay-1) - play(icol,ilay)));
        if (ilay == nlay+1) {
          tlev_wk(icol,nlay+1) = tlay(icol,nlay) + (plev(icol,nlay+1)-play(icol,nlay))*(tlay(icol,nlay)-tlay(icol,nlay-1)) / 
                                                   (play(icol,nlay)-play(icol,nlay-1));
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
    parallel_for( Bounds<2>(ngpt,ncol) , YAKL_LAMBDA (int igpt, int icol) {
      sources_sfc_source(icol,igpt) = sfc_source_t(igpt,icol);
    });
    reorder123x321(ngpt, nlay, ncol, lay_source_t    , sources.lay_source    );
    reorder123x321(ngpt, nlay, ncol, lev_source_inc_t, sources.lev_source_inc);
    reorder123x321(ngpt, nlay, ncol, lev_source_dec_t, sources.lev_source_dec);
  }



  // Utility function, provided for user convenience
  // computes column amounts of dry air using hydrostatic equation
  real2d get_col_dry(real2d const &vmr_h2o, real2d const &plev, real1d const &latitude=real1d()) {
    // first and second term of Helmert formula
    real constexpr helmert1 = 9.80665_wp;
    real constexpr helmert2 = 0.02586_wp;
    int ncol = size(plev,1);
    int nlev = size(plev,2);
    real1d g0("g0",size(plev,1));
    if (allocated(latitude)) {
      // A purely OpenACC implementation would probably compute g0 within the kernel below
      // do icol = 1, ncol
      parallel_for( Bounds<1>(ncol) , YAKL_LAMBDA (int icol) {
        g0(icol) = helmert1 - helmert2 * cos(2.0_wp * M_PI * latitude(icol) / 180.0_wp); // acceleration due to gravity [m/s^2]
      });
    } else {
      // do icol = 1, ncol
      auto &grav = ::grav;
      parallel_for( Bounds<1>(ncol) , YAKL_LAMBDA (int icol) {
        g0(icol) = grav;
      });
    }

    real2d col_dry("col_dry",size(plev,1),size(plev,2)-1);
    // do ilev = 1, nlev-1
    //   do icol = 1, ncol
    auto &m_dry = ::m_dry;
    parallel_for( Bounds<2>(nlev-1,ncol) , YAKL_LAMBDA (int ilev , int icol) {
      real delta_plev = abs(plev(icol,ilev) - plev(icol,ilev+1));
      // Get average mass of moist air per mole of moist air
      real fact = 1._wp / (1.+vmr_h2o(icol,ilev));
      real m_air = (m_dry + m_h2o * vmr_h2o(icol,ilev)) * fact;
      col_dry(icol,ilev) = 10._wp * delta_plev * avogad * fact/(1000._wp*m_air*100._wp*g0(icol));
    });
    return col_dry;
  }



 // Utility function to combine optical depths from gas absorption and Rayleigh scattering
 //   (and reorder them for convenience, while we're at it)
 void combine_and_reorder(real3d const &tau, real3d const &tau_rayleigh, bool has_rayleigh, OpticalProps1scl &optical_props) {
    int ncol = size(tau,3);
    int nlay = size(tau,2);
    int ngpt = size(tau,1);
    reorder123x321(ngpt, nlay, ncol, tau, optical_props.tau);
  }



 // Utility function to combine optical depths from gas absorption and Rayleigh scattering
 //   (and reorder them for convenience, while we're at it)
 void combine_and_reorder(real3d const &tau, real3d const &tau_rayleigh, bool has_rayleigh, OpticalProps2str &optical_props) {
    int ncol = size(tau,3);
    int nlay = size(tau,2);
    int ngpt = size(tau,1);
    if (has_rayleigh) {
      // combine optical depth and rayleigh scattering
      combine_and_reorder_2str(ncol, nlay, ngpt, tau, tau_rayleigh, optical_props.tau, optical_props.ssa, optical_props.g);
    } else {
      // index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
      reorder123x321(ngpt, nlay, ncol, tau, optical_props.tau);
      zero_array(optical_props.ssa);
      zero_array(optical_props.g  );
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

