
#pragma once

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
// Encapsulate source function arrays for longwave/lw/internal sources
//    and shortwave/sw/external source.
//
// -------------------------------------------------------------------------------------------------

// Type for longwave sources: computed at layer center, at layer edges using
//   spectral mapping in each direction separately, and at the surface
// Not implementing get_subset because it isn't used
class SourceFuncLW : public OpticalProps {
public:
  real3d lay_source;
  real3d lev_source_inc;
  real3d lev_source_dec;
  real2d sfc_source;


  bool is_allocated() const  { return this->is_initialized() && yakl::intrinsics::allocated(this->sfc_source); }


  void alloc(int ncol, int nlay) {
    if (! this->is_initialized()) { stoprun("source_func_lw%alloc: not initialized so can't allocate"); }
    if (ncol <= 0 || nlay <= 0) { stoprun("source_func_lw%alloc: must provide positive extents for ncol, nlay"); }
    int ngpt = this->get_ngpt();
    this->sfc_source     = real2d("sfc_source"    ,ncol,ngpt);
    this->lay_source     = real3d("lay_source"    ,ncol,nlay,ngpt);
    this->lev_source_inc = real3d("lev_source_inc",ncol,nlay,ngpt);
    this->lev_source_dec = real3d("lev_source_dec",ncol,nlay,ngpt);
  }


  void alloc(int ncol, int nlay, OpticalProps const &op) {
    if (! op.is_initialized()) { stoprun("source_func_lw::alloc: op not initialized"); }
    this->finalize();
    this->init(op);
    this->alloc(ncol,nlay);
  }


  void finalize() {
    this->lay_source     = real3d();
    this->lev_source_inc = real3d();
    this->lev_source_dec = real3d();
    this->sfc_source     = real2d();
    OpticalProps::finalize();
  }


  int get_ncol() const {
    if (this->is_allocated()) { return yakl::intrinsics::size(this->lay_source,1); } else { return 0; }
  }


  int get_nlay() const {
    if (this->is_allocated()) { return yakl::intrinsics::size(this->lay_source,2); } else { return 0; }
  }


  void print_norms() const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

                                     std::cout << "name          : " << name                << "\n";
    if (allocated(lay_source    )) { std::cout << "lay_source    : " << sum(lay_source    ) << "\n"; }
    if (allocated(lev_source_inc)) { std::cout << "lev_source_inc: " << sum(lev_source_inc) << "\n"; }
    if (allocated(lev_source_dec)) { std::cout << "lev_source_dec: " << sum(lev_source_dec) << "\n"; }
    if (allocated(sfc_source    )) { std::cout << "sfc_source    : " << sum(sfc_source    ) << "\n"; }
    if (allocated(band2gpt      )) { std::cout << "band2gpt      : " << sum(band2gpt      ) << "\n"; }
    if (allocated(gpt2band      )) { std::cout << "gpt2band      : " << sum(gpt2band      ) << "\n"; }
    if (allocated(band_lims_wvn )) { std::cout << "band_lims_wvn : " << sum(band_lims_wvn ) << "\n"; }
  }

};



// Type for shortave sources: top-of-domain spectrally-resolved flux
// Not implementing get_subset because it isn't used
class SourceFuncSW : public OpticalProps {
public:
  real2d toa_source;


  bool is_allocated() const { return this->is_initialized() && yakl::intrinsics::allocated(this->toa_source); }


  void alloc(int ncol) {
    if (! this->is_initialized()) { stoprun("source_func_sw%alloc: not initialized so can't allocate"); }
    if (ncol <= 0) { stoprun("source_func_sw%alloc: must provide positive extents for ncol"); }
    this->toa_source = real2d("toa_source",ncol,this->get_ngpt());
  }


  void alloc(int ncol, OpticalProps const &op) {
    if (! op.is_initialized()) { stoprun("source_func_sw::alloc: op not initialized"); }
    this->init(op);
    this->alloc(ncol);
  }


  void finalize() {
    this->toa_source = real2d();
    OpticalProps::finalize();
  }


  int get_ncol() const {
    if (this->is_allocated()) { return yakl::intrinsics::size(this->toa_source,1); } else { return 0; }
  }


  void print_norms() const {
    using yakl::intrinsics::sum;
    using yakl::intrinsics::allocated;

                                     std::cout << "name          : " << name                << "\n";
    if (allocated(toa_source    )) { std::cout << "toa_source    : " << sum(toa_source    ) << "\n"; }
    if (allocated(band2gpt      )) { std::cout << "band2gpt      : " << sum(band2gpt      ) << "\n"; }
    if (allocated(gpt2band      )) { std::cout << "gpt2band      : " << sum(gpt2band      ) << "\n"; }
    if (allocated(band_lims_wvn )) { std::cout << "band_lims_wvn : " << sum(band_lims_wvn ) << "\n"; }
  }

};

