#pragma once

#include "rrtmgp_const.h"
#include "mo_fluxes_broadband_kernels.h"
#include <iomanip>

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
// Compute output quantities from RTE based on spectrally-resolved flux profiles
//    This module contains an abstract class and a broadband implmentation that sums over all spectral points
//    The abstract base class defines the routines that extenstions must implement: reduce() and are_desired()
//    The intent is for users to extend it as required, using mo_flxues_broadband as an example
//
// -------------------------------------------------------------------------------------------------

class FluxesBroadband {
public:
  real2d flux_up;
  real2d flux_dn;
  real2d flux_net;
  real2d flux_dn_dir;


  void reduce(real3d const &gpt_flux_up, const real3d &gpt_flux_dn, OpticalProps const &spectral_disc,
              bool top_at_1, real3d const &gpt_flux_dn_dir=real3d()) {
    using yakl::intrinsics::size;
    using yakl::intrinsics::allocated;

    int ncol = size(gpt_flux_up,1);
    int nlev = size(gpt_flux_up,2);
    int ngpt = size(gpt_flux_up,3);

    // Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    if (allocated(this->flux_dn_dir) && ! allocated(gpt_flux_dn_dir)) {
      stoprun("reduce: requesting direct downward flux but this hasn't been supplied");
    }

    // Broadband fluxes - call the kernels
    if (allocated(this->flux_up    )) { sum_broadband(ncol, nlev, ngpt, gpt_flux_up,     this->flux_up    ); }
    if (allocated(this->flux_dn    )) { sum_broadband(ncol, nlev, ngpt, gpt_flux_dn,     this->flux_dn    ); }
    if (allocated(this->flux_dn_dir)) { sum_broadband(ncol, nlev, ngpt, gpt_flux_dn_dir, this->flux_dn_dir); }
    if (allocated(this->flux_net   )) {
      // Reuse down and up results if possible
      if (allocated(this->flux_dn) && allocated(this->flux_up)) {
        net_broadband(ncol, nlev,      this->flux_dn, this->flux_up, this->flux_net);
      } else {
        net_broadband(ncol, nlev, ngpt,  gpt_flux_dn,   gpt_flux_up, this->flux_net);
      }
    }
  }


  bool are_desired() const {
    using yakl::intrinsics::allocated;
    return allocated(this->flux_up) || allocated(this->flux_dn) || allocated(this->flux_dn_dir) || allocated(this->flux_net);
  }


  void print_norms() const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

    if (allocated(flux_up    )) { std::cout << std::setprecision(16) << "flux_up    : " << sum(flux_up    ) << "\n"; }
    if (allocated(flux_dn    )) { std::cout << std::setprecision(16) << "flux_dn    : " << sum(flux_dn    ) << "\n"; }
    if (allocated(flux_net   )) { std::cout << std::setprecision(16) << "flux_net   : " << sum(flux_net   ) << "\n"; }
    if (allocated(flux_dn_dir)) { std::cout << std::setprecision(16) << "flux_dn_dir: " << sum(flux_dn_dir) << "\n"; }
  }

};

#ifdef RRTMGP_ENABLE_KOKKOS
class FluxesBroadbandK {
public:
  real2dk flux_up;
  real2dk flux_dn;
  real2dk flux_net;
  real2dk flux_dn_dir;


  void reduce(real3dk const &gpt_flux_up, const real3dk &gpt_flux_dn, OpticalPropsK const &spectral_disc,
              bool top_at_1, real3dk const &gpt_flux_dn_dir=real3dk()) {
    int ncol = gpt_flux_up.extent(0);
    int nlev = gpt_flux_up.extent(1);
    int ngpt = gpt_flux_up.extent(2);

    // Self-consistency -- shouldn't be asking for direct beam flux if it isn't supplied
    if (this->flux_dn_dir.is_allocated() && ! gpt_flux_dn_dir.is_allocated()) {
      stoprun("reduce: requesting direct downward flux but this hasn't been supplied");
    }

    // Broadband fluxes - call the kernels
    if (this->flux_up.is_allocated()    ) { sum_broadband(ncol, nlev, ngpt, gpt_flux_up,     this->flux_up    ); }
    if (this->flux_dn.is_allocated()    ) { sum_broadband(ncol, nlev, ngpt, gpt_flux_dn,     this->flux_dn    ); }
    if (this->flux_dn_dir.is_allocated()) { sum_broadband(ncol, nlev, ngpt, gpt_flux_dn_dir, this->flux_dn_dir); }
    if (this->flux_net.is_allocated()   ) {
      // Reuse down and up results if possible
      if (this->flux_dn.is_allocated() && this->flux_up.is_allocated()) {
        net_broadband(ncol, nlev,      this->flux_dn, this->flux_up, this->flux_net);
      } else {
        net_broadband(ncol, nlev, ngpt,  gpt_flux_dn,   gpt_flux_up, this->flux_net);
      }
    }
  }


  bool are_desired() const {
    return this->flux_up.is_allocated() || this->flux_dn.is_allocated() || this->flux_dn_dir.is_allocated() || this->flux_net.is_allocated();
  }


  void print_norms() const {
    if (flux_up.is_allocated()    ) { std::cout << std::setprecision(16) << "flux_up    : " << conv::sum(flux_up    ) << "\n"; }
    if (flux_dn.is_allocated()    ) { std::cout << std::setprecision(16) << "flux_dn    : " << conv::sum(flux_dn    ) << "\n"; }
    if (flux_net.is_allocated()   ) { std::cout << std::setprecision(16) << "flux_net   : " << conv::sum(flux_net   ) << "\n"; }
    if (flux_dn_dir.is_allocated()) { std::cout << std::setprecision(16) << "flux_dn_dir: " << conv::sum(flux_dn_dir) << "\n"; }
  }

  void validate_kokkos(const FluxesBroadband& orig)
  {
    conv::compare_yakl_to_kokkos(orig.flux_up, flux_up);
    conv::compare_yakl_to_kokkos(orig.flux_dn, flux_dn);
    conv::compare_yakl_to_kokkos(orig.flux_net, flux_net);
    conv::compare_yakl_to_kokkos(orig.flux_dn_dir, flux_dn_dir);
  }

};
#endif
